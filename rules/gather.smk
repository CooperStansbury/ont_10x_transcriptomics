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
        refgenome=config['reference']['fasta'],
    output:
        OUTPUT_PATH + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"


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



rule get_annotations:
    """
    Copies the annotation file to the references directory.
    """
    input:
        config['reference']['annotation'],
    output:
        OUTPUT_PATH + 'references/annotations.gtf'
    shell:
        "cp {input} {output}"

