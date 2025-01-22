rule get_reference:
    input:
        refgenome=config['reference']['fasta'],
    output:
        OUTPUT_PATH + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"


rule get_reference_genes:
    input:
        config['transcript_path'],
    output:
        OUTPUT_PATH + 'references/transcripts.fa.gz'
    shell:
        "cp {input} {output}"


rule get_annotations:
    input:
        refgenome=config['gtf_path'],
    output:
        OUTPUT_PATH + 'references/annotations.gtf'
    shell:
        "cp {input} {output}"


        
            
rule get_gene_table:
    input:
        annotations=config['gtf_path'],
    output:
        OUTPUT_PATH + "references/geneTable.csv"
    conda:
        "bioinf"
    shell:
        "python scripts/getGeneTable.py {input} {output}"


rule minimap2_index:
    input:
        refgenome=OUTPUT_PATH + 'references/reference.fa.gz'
    output:
        OUTPUT_PATH + 'references/reference.mmi'
    threads:
        config['threads']
    conda:
        "aligner"
    shell:
        "minimap2 -t {threads} -d {output} {input.refgenome}"

