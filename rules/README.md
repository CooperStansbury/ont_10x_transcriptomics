# Snakemake Workflow

This directory contains Snakemake rule files (`.smk`) that define the steps for processing genomic data.

## Workflow Stages

The workflow is divided into two main stages:

### 1. Data Preparation (`gather`)

This stage organizes and prepares your data and reference files:

*   Copies configuration files, input lists, reference genome, annotations, and barcode whitelists.
*   Organizes files into the correct project directories.
*   Extracts genome build information.

### 2. Mapping and Analysis (`core`)

This stage processes your data:

*   Indexes the reference genome using `minimap2`.
*   Generates quality control reports for raw and demultiplexed reads using `seqkit`.
*   Demultiplexes FASTQ files based on cell barcodes using `blaze`.
*   Aligns demultiplexed reads to the reference genome using `minimap2`.
*   Filters alignments based on mapping quality using `samtools`.
*   Tags BAM files with barcode, UMI, and read name information.
*   Aggregates reads by chromosome and sorts them using `samtools`.
*   Generates chromosome specific BAM files
*   Counts reads overlapping genomic features (e.g., genes) using `HTSeq-count` for each chromosome.
*   Compiles and annotates a single AnnData object from multiple chromosome-specific h5ad files, including gene filtering and sparse matrix conversion.

### 3. Reporting (`reporting`)

This stage generates quality control and summary reports:

*   Runs `NanoQC` on FASTQ files to generate quality control reports.
*   Generates sequencing statistics from FASTQ files using `NanoStat`.
*   Creates summary reports of key statistics for raw and demultiplexed FASTQ files using `seqkit`.
*   Generates alignment summary statistics using `samtools flagstat`.

## Output

The workflow produces:

*   Indexed reference genome (`.mmi`).
*   Raw FASTQ report (`raw_fastq_report.txt`).
*   Demultiplexed FASTQ report (`demultiplexed_fastq_report.txt`).
*   Demultiplexing summary files, including matched reads, empty barcode lists, knee plots, putative barcode lists, and summary reports.
*   Aligned and indexed BAM files for each sample (`.bam`, `.bam.bai`).
*   Chromosome-specific, sorted, and indexed BAM files.
*   Tagged BAM files with barcode and UMI information (`.tagged.bam`).
*   CSV files containing alignment records with barcode and UMI information (`.records.csv`).
*   Chromosome-specific read counts in h5ad format (`.counts.h5ad`).
*   A single, compiled, and annotated AnnData object (`anndata.raw.h5ad`).
*   NanoQC reports for each sample.
*   NanoStat reports for each sample.
*   Alignment summary statistics for each sample (`.flagstat.txt`).
