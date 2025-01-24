# This Snakemake workflow manages indexing a reference genome and annotations,
# prepares input FASTQ files, and sets up configuration files.


rule get_config:
    """
    Copies the configuration file to the output directory.

    This rule takes the configuration file specified by `CONFIG_ABS_PATH`
    and copies it to the `config/` directory in the output location.
    """
    input:
        CONFIG_ABS_PATH
    output:
        OUTPUT_PATH + "config/" + CONFIG_BASENAME,
    shell:
        """ cp {input} {output} """


rule get_input_list:
    """
    Copies the input file list to the output directory.

    This rule takes the input file list specified by `INPUT_ABS_PATH`
    and copies it to the `config/` directory in the output location. 
    This list likely contains information about the FASTQ files 
    to be processed.
    """
    input:
        INPUT_ABS_PATH
    output:
        OUTPUT_PATH + "config/" + INPUT_BASENAME,
    shell:
        """ cp {input} {output} """
    


rule get_fastq_files:
    """
    Copies input FASTQ files to the output directory.

    This rule takes a list of input FASTQ file paths and copies them 
    to the specified output paths. This is typically used to stage 
    input data for downstream processing.
    """
    input:
        input_file_paths
    output:
        output_file_paths
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input):
            outPath = output[i]
            copyfile(refPath, outPath)


rule get_reference:
    """
    Copies the reference genome to the references directory.
    """
    input:
        refgenome=config['inputs']['reference_fasta'],
    output:
        OUTPUT_PATH + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"


rule get_whitelist:
    """
    Copy the barcode whitelist file to the references directory.
    """
    input:
        config['demux']['barcode_whitelist']
    output:
        OUTPUT_PATH + "references/barcode_whitelist.txt",
    shell:
        """ cp {input} {output} """


rule get_annotations:
    """
    Copies the annotation file to the references directory.
    """
    input:
        config['inputs']['reference_annotation'],
    output:
        OUTPUT_PATH + 'references/annotations.gtf'
    shell:
        "cp {input} {output}"


rule get_gene_table:
    """
    Extracts gene information from a GTF file and creates a gene table,
    excluding entries without a gene_name.

    The output table contains the following columns:
    - gene_id
    - gene_name
    - gene_biotype
    - chrom
    - start
    - end
    """
    input:
        OUTPUT_PATH + 'references/annotations.gtf',
    output:
        OUTPUT_PATH + 'references/gene_table.tsv',
    shell:
        """
        awk '
            BEGIN {{ FS = "\t" }}
            {{
                if ($3 == "gene") {{
                    split($9, a, ";");
                    gene_id = "";
                    gene_name = "";
                    gene_biotype = "";

                    for (i in a) {{
                        if (a[i] ~ /gene_id/) {{
                            split(a[i], b, " ");
                            gene_id = b[2];
                        }}
                        if (a[i] ~ /gene_name/) {{
                            split(a[i], b, " ");
                            gene_name = b[2];
                        }}
                        if (a[i] ~ /gene_biotype/) {{
                            split(a[i], b, " ");
                            gene_biotype = b[2];
                        }}
                    }}

                    gsub(/"/, "", gene_id);
                    gsub(/"/, "", gene_name);
                    gsub(/"/, "", gene_biotype);

                    if (gene_name != "") {{
                        print gene_id "\t" gene_name "\t" gene_biotype "\t" $1 "\t" $4 "\t" $5
                    }}
                }}
            }}
        ' {input} | sed '1igene_id\tgene_name\tgene_biotype\tchrom\tstart\tend' > {output}
        """


rule get_genome_build:
    """
    Extracts all comment lines (starting with '#!') from a reference GTF file.

    This rule uses grep to find all lines beginning with '#!' within 
    the GTF file and saves them to a separate output file.
    """
    input:
        gtf=OUTPUT_PATH + 'references/annotations.gtf' 
    output:
        build=OUTPUT_PATH + 'references/genome_build.txt'
    shell:
        """
        grep "^#!" {input} > {output}
        """


rule get_chroms:
    """
    Extracts unique chromosome names from a GTF file, excluding those with ".".

    This is a hacky way to ignore unplaced contigs.

    This rule uses command-line tools to identify and extract 
    the unique chromosome names present in the input GTF file,
    filtering out any chromosome names containing a period.
    """
    input:
        gtf=OUTPUT_PATH + 'references/annotations.gtf'
    output:
        chroms=OUTPUT_PATH + 'references/chroms.txt'
    shell:
        """
        cut -f 1 {input.gtf} | grep -v '^#' | sort -u | grep -v '\.' > {output.chroms}
        """

""" DEFINE THE CHROMOSOMES """
chromosomes = pu.get_names(rules.get_chroms.output[0])

rule create_chromosome_gtf:
    """
    This rule extracts GTF entries for a specific chromosome from the main
    annotations.gtf file and creates a chromosome-specific GTF file.
    """
    input:
        gtf=OUTPUT_PATH + "references/annotations.gtf",
    params:
        chrom=lambda wildcards: wildcards.chrom
    output:
        OUTPUT_PATH + "references/by_chrom/{chrom}.gtf"
    shell:
        """
        grep "^{params.chrom}\s" {input.gtf} > {output}
        """