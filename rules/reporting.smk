rule nanoqc:
    """
    Runs NanoQC on a FASTQ file to generate quality control reports.
    """
    input:
        fastq=OUTPUT_PATH + "fastq/{sid}.fastq" + extension
    output:
        OUTPUT_PATH + "reports/nanoqc/{sid}/{sid}.done"
    params:
        output_dir=OUTPUT_PATH + "reports/nanoqc/{sid}"
    conda:
        "../envs/pipeline-QC.yaml"
    shell:
        """
        nanoQC -o {params.output_dir} {input.fastq}
        touch {output}  
        """


rule nanostat:
    """
    Generates sequencing statistics from a FASTQ file using NanoStat.
    """
    input:
        fastq=OUTPUT_PATH + "fastq/{sid}.fastq" + extension
    output:
        OUTPUT_PATH + "reports/nanostat/{sid}.txt"
    conda:
        "../envs/pipeline-QC.yaml"
    shell:
        """
        NanoStat --fastq {input.fastq} > {output}
        """


rule raw_report:
    """
    Generates a summary report of key statistics for raw FASTQ files using seqkit.
    """
    input:
        expand(OUTPUT_PATH + "fastq/{sid}.fastq" + extension, sid=samples),
    output:
        OUTPUT_PATH + "reports/seqkit_stats/raw_fastq_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "../envs/pipeline-QC.yaml"
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""



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
        "../envs/pipeline-QC.yaml"
    params:
        files=expand(OUTPUT_PATH + "demultiplex/{sid}.matched_reads.fastq" + extension, sid=samples),
    shell:
        """seqkit stats -a -b -j {threads} {params.files} -o {output}"""


rule alignment_summary:
    """Generates alignment summary statistics using samtools flagstat.
    
    This rule takes a tagged BAM file as input and produces a text file
    containing alignment statistics, such as the number of mapped and unmapped reads,
    duplicate reads, and properly paired reads (for paired-end data).
    """
    input:
        OUTPUT_PATH + 'mapping/{sid}.tagged.bam',
    output:
         OUTPUT_PATH + "reports/alignment/{sid}.flagstat.txt",
    conda:
        "../envs/pipeline-QC.yaml"
    shell:
        """samtools flagstat {input} > {output}"""
