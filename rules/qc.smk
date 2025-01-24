rule longqc_ont:
    """
    Performs quality control (QC) on ONT sequencing data using LongQC.
    Runs LongQC in ONT mode to generate quality control statistics 
    for each sample’s raw FASTQ file.
    """
    input:
        fastq=OUTPUT_PATH + "fastq/{sid}" + extension
    output:
        report_dir=directory(OUTPUT_PATH + "qc/longqc/{sid}/"),
        summary=OUTPUT_PATH + "qc/longqc/{sid}/summary.txt"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)])
    threads:
        config['threads'] // 4
    conda:
        "pipeline-core"
    shell:
        """
        longQC.py sampleqc -x ont -s {wildcards.sid} -o {output.report_dir} {input.fastq}
        mv {output.report_dir}/summary.txt {output.summary}
        """