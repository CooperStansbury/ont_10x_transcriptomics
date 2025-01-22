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