# Snakemake Workflow

This directory contains Snakemake rule files (`.smk`) that define the steps for processing genomic data.

## Workflow Stages

The workflow is divided into two main stages:

### 1. Data Preparation (`gather`)

This stage organizes and prepares your data and reference files:

* Copies configuration files, input lists, reference genome, annotations, and barcode whitelists.
* Organizes files into the correct project directories.
* Extracts genome build information.

### 2. Mapping and Analysis (`core`)

This stage processes your data:

* Indexes the reference genome.
* Generates quality control reports for raw and processed reads.
* Demultiplexes FASTQ files.
* Aligns reads to the reference genome.
* Tags and indexes BAM files.
* Generates alignment statistics.

## Output

The workflow produces:

* Aligned and indexed BAM files.
* Summary reports with quality control and alignment statistics. 