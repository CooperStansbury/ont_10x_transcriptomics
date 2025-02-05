rule minimap2_index:
    """
    Indexes the reference genome using minimap2.
    """
    input:
        OUTPUT_PATH + 'references/reference.fa.gz'
    output:
        OUTPUT_PATH + 'references/reference.mmi'
    threads:
        config['threads']
    conda:
        "../envs/pipeline-core.yaml"
    shell:
        "minimap2 -t {threads} -d {output} {input}"


rule demultiplex:
    """
    Demultiplexes FASTQ files based on cell barcodes.
    """
    input:
        fastq=OUTPUT_PATH + "fastq/{sid}.fastq" + extension,
        whitelist=config['demux']['barcode_whitelist'],
    output:
        OUTPUT_PATH + 'demultiplex/{sid}.matched_reads.fastq' + extension,
        OUTPUT_PATH + 'demultiplex/{sid}.emtpy_bc_list.csv',
        OUTPUT_PATH + 'demultiplex/{sid}.knee_plot.png',
        OUTPUT_PATH + 'demultiplex/{sid}.putative_bc.csv',
        OUTPUT_PATH + 'demultiplex/{sid}.summary.txt',
        OUTPUT_PATH + 'demultiplex/{sid}.whitelist.csv',
        touch(OUTPUT_PATH + "demultiplex/{sid}.done"),
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] 
    params:
        expected=config['demux']['expected_cells'],
        output_prefix=lambda wildcards: OUTPUT_PATH + "demultiplex/" + wildcards.sid + ".", 
    conda:
        "../envs/pipeline-core.yaml"
    log:
        OUTPUT_PATH + "logs/demultiplex/{sid}.log",
    shell:
        """blaze --expect-cells {params.expected} \
        --output-prefix {params.output_prefix} --threads {threads} \
        --full-bc-whitelist {input.whitelist} {input.fastq} > {log}"""


rule align_reads:
    """
    Align demultiplexed reads to a reference genome using minimap2, filter by mapping quality, and sort the output using samtools.

    This rule takes demultiplexed reads and aligns them to a reference genome using minimap2 with specified parameters. 
    The resulting SAM output is then filtered based on the provided mapping quality threshold using samtools view.
    Finally, the filtered SAM output is sorted and converted to BAM format using samtools.
    """
    input:
        flag=OUTPUT_PATH + "demultiplex/{sid}.done",
        ref=OUTPUT_PATH + 'references/reference.fa.gz',
        refindex=OUTPUT_PATH + 'references/reference.mmi',
    output:
        bam=OUTPUT_PATH + 'mapping/{sid}.bam',
    params:
        args=config['alignment']['minimap2_params'],
        mapq=config['alignment']['mapq_threshold'],
        fastq=OUTPUT_PATH + "demultiplex/{sid}.matched_reads.fastq" + extension,
    threads:
        int(config['threads'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    log:
        OUTPUT_PATH + "logs/mapping/{sid}.log"
    conda:
        "../envs/pipeline-core.yaml"
    shell:
        """
        minimap2 {params.args} -t {threads} \
        {input.ref} {params.fastq} | \
        samtools view -h -q {params.mapq} - | \
        samtools sort -@ {threads} -O bam -o {output.bam} 2>&1 | tee {log}
        """


rule tag_bam:
    """
    Tag BAM files with barcode, UMI, and read name information, and save relevant information in a CSV file.

    This rule takes a BAM file, tags each alignment with the barcode, UMI, and read name extracted from the query name,
    and writes the modified alignments to a new BAM file. It also stores relevant information about each alignment in a CSV file.
    """
    input:
        OUTPUT_PATH + 'mapping/{sid}.bam'
    output:
        bam=OUTPUT_PATH + 'mapping/{sid}.tagged.bam',
        records=OUTPUT_PATH + 'mapping/{sid}.records.csv',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "../envs/pipeline-core.yaml"
    shell:
        """python scripts/tag_bam.py {input} {output.bam} {output.records}"""



rule aggregate_reads_by_chromosome:
    """
    Aggregates reads from multiple BAM files into a single BAM file for each chromosome.

    This rule takes a list of BAM files (typically one per sample) and merges them.
    Then, it filters and sorts the merged reads by chromosome, producing one
    chromosome-specific, sorted BAM file, and its corresponding index.
    """
    input:
        bam_files=expand(OUTPUT_PATH + 'mapping/{sid}.tagged.bam', sid=samples),
    params:
        chrom=lambda wc: wc.chrom
    output:
        bam=OUTPUT_PATH + "mapping/by_chrom/{chrom}.bam",
        bai=OUTPUT_PATH + "mapping/by_chrom/{chrom}.bam.bai"
    conda:
        "../envs/pipeline-core.yaml"
    threads:
        int(config['threads'])
    shell:
        """
        samtools merge -@ {threads} - {input.bam_files} | \
        samtools view -b -@ {threads} -h -r {params.chrom} - | \
        samtools sort -@ {threads} - -o {output.bam}
        samtools index {output.bam} {output.bai}
        """


rule htseq_count:
    """
    Counts reads overlapping genomic features (e.g., genes) using HTSeq-count,
    processing each chromosome separately.

    This rule takes a chromosome-specific BAM file and a corresponding GTF
    file as input and produces a counts file using HTSeq-count with barcode
    correction. The counts file lists the number of reads mapping to each
    feature in the GTF file. This rule processes each chromosome independently
    for improved parallelism.
    """
    input:
        bam=OUTPUT_PATH + "mapping/by_chrom/{chrom}.bam",
        annotations=OUTPUT_PATH + "references/by_chrom/{chrom}.gtf"
    output:
        OUTPUT_PATH + "counts/{chrom}.counts.h5ad"
    conda:
        "../envs/pipeline-core.yaml"
    params:
        d=lambda wc: int(config['counts']['umi_distance']),
        chrom=lambda wc: wc.chrom
    shell:
        """
        htseq-count-barcodes \
            --correct-UMI-distance {params.d} \
            --counts_output_sparse \
            {input.bam} \
            {input.annotations} -c {output}
        """


def make_counts_output(wildcards):
    """
    Generates a list of chromosome-specific GTF file paths based on the output of the get_chroms checkpoint.
    """
    chromosomes = [line.strip() for line in open(checkpoints.get_chroms.get().output.chroms)]
    return expand(OUTPUT_PATH + "counts/{chrom}.counts.h5ad", chrom=chromosomes)


rule compile_anndata:
    """
    Compiles and annotates a single AnnData object from multiple chromosome-specific h5ad files.

    This rule takes a list of h5ad files (each typically representing counts for a single chromosome)
    and a gene metadata table as input. It uses the 'compile_anndata.py' script to perform the following:

    1. Concatenation: Combines the individual h5ad files into a single AnnData object.
    2. Annotation: Merges the provided gene metadata table with the AnnData object's `var` attribute.
    3. Gene Filtering: Removes genes that lack a 'gene_name' entry in the metadata.
    4. Sparse Matrix Conversion: Ensures the resulting count matrix (`adata.X`) is stored as a sparse 
       matrix (CSR format) to optimize memory usage.

    The output of this rule is a single, comprehensive h5ad file containing the combined and annotated
    count data, ready for downstream analysis.
    """
    input:
        h5ad=make_counts_output,
        gene_table=OUTPUT_PATH + 'references/gene_table.tsv',
    output:
        OUTPUT_PATH + 'anndata/anndata.raw.h5ad',
    conda:
        "../envs/pipeline-core.yaml"
    shell:
        """ python scripts/compile_anndata.py {output} {input.gene_table} {input.h5ad}"""


rule generate_whitelist:
    """
    Generates a whitelist of barcodes from a Scanpy object.
    """
    input:
        OUTPUT_PATH + 'anndata/anndata.raw.h5ad',
    output:
        OUTPUT_PATH + 'whitelist/detected_barcodes.txt'
    conda:
        "../envs/pipeline-core.yaml"
    shell:
        """ python scripts/get_whitelist.py {input} {output} """
        
    
