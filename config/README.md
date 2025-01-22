# Snakemake Pipeline Configuration

This directory contains configuration files used to set up and run a Snakemake pipeline. The configuration files specify parameters, input/output paths, and execution settings required for different stages of the workflow.

## Overview

- `config.yaml` – Main configuration file with key parameters for the pipeline.


## Config

This YAML configuration file provides settings for running the pipeline, including resource allocation, input files, and reference genome information.

### Parameters

- **Threads:** `36`  
  Number of threads to be used for processing.

- **Output Path:**  
  `/scratch/indikar_root/indikar1/shared_data/ont_10x_pipeline_test/`  
  Directory where pipeline outputs will be stored.

- **FASTQ Paths:**  
  `config/fastq_paths.txt`  
  Path to the text file containing a list of input FASTQ files.

- **FASTA File:**  
  `/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.cdna.all.fa`  
  Reference genome sequence in FASTA format. Expects UNZIPPED `.fa` file.

- **Annotation File:**  
  `/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.107.gtf`  
  Gene annotation file in GTF format. Expects UNZIPPED `.gtf` file.



