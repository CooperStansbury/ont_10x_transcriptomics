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
        "pipeline-core"
    shell:
        "minimap2 -t {threads} -d {output} {input}"


rule raw_report:
    """
    Generates a summary report of key statistics for raw FASTQ files using seqkit.
    """
    input:
        expand(OUTPUT_PATH + "fastq/{sid}" + extension, sid=samples),
    output:
        OUTPUT_PATH + "reports/seqkit_stats/raw_fastq_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "pipeline-core"
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule demultiplex:
    """
    Demultiplexes FASTQ files based on cell barcodes.
    """
    input:
        fastq=OUTPUT_PATH + "fastq/{sid}" + extension,
        whitelist=config['demux']['barcode_whitelist'],
    output:
        OUTPUT_PATH + 'demultiplex/{sid}.matched_reads' + extension,
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
        "pipeline-core"
    log:
        OUTPUT_PATH + "logs/demultiplex/{sid}.log",
    shell:
        """blaze --expect-cells {params.expected} \
        --output-prefix {params.output_prefix} --threads {threads} \
        --full-bc-whitelist {input.whitelist} {input.fastq} > {log}"""


rule demultiplexed_report:
    """
    Generates a summary report of key statistics for demultiplexed FASTQ files using seqkit.
    """
    input:
        flags=expand(OUTPUT_PATH + "demultiplex/{sid}.done", sid=samples),
    output:
        OUTPUT_PATH + "reports/seqkit_stats/demultiplexed_fastq_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "pipeline-core"
    params:
        files=expand(OUTPUT_PATH + "demultiplex/{sid}.matched_reads" + extension, sid=samples),
    shell:
        """seqkit stats -a -b -j {threads} {params.files} -o {output}"""



rule align_reads:
    """
    Align demultiplexed reads to a reference genome using minimap2 and sort the output using samtools.

    This rule takes demultiplexed reads and aligns them to a reference genome using minimap2 with specified parameters. 
    The resulting SAM output is then sorted and converted to BAM format using samtools.
    """
    input:
        flag=OUTPUT_PATH + "demultiplex/{sid}.done",
        ref=OUTPUT_PATH + 'references/reference.fa.gz',
        refindex=OUTPUT_PATH + 'references/reference.mmi',
    output:        
        bam=OUTPUT_PATH + 'mapping/{sid}.bam',
    params:
        args=config['minimap2_params'],
        fastq=OUTPUT_PATH + "demultiplex/{sid}.matched_reads" + extension,
    threads:
        int(config['threads'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    log:
        OUTPUT_PATH + "logs/mapping/{sid}.log"
    conda:
        "pipeline-core"
    shell:
        """minimap2 {params.args} -t {threads} \
        {input.ref} {params.fastq} | samtools sort \
        -@ {threads} -O bam -o {output.bam} 2>&1 | tee {log}"""


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
        "pipeline-core"
    shell:
        """python scripts/tag_bam.py {input} {output.bam} {output.records}"""


rule samtools_index:
    """
    Index a BAM file using samtools.

    This rule creates a BAM index (.bai) file for a given BAM file using samtools. 
    The index allows for efficient random access to the BAM file.
    """
    input:
        OUTPUT_PATH + 'mapping/{sid}.tagged.bam'
    output:
        OUTPUT_PATH + 'mapping/{sid}.tagged.bam.bai'
    conda:
        "pipeline-core"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    shell:
        """samtools index -@ {threads} {input}"""

