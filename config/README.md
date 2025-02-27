# Pipeline Configuration

This directory contains configuration files for setting up and running a Snakemake pipeline. These files define essential parameters, input/output paths, and execution settings for different stages of the workflow.

## Overview

The following configuration files are included:

-   **`config.yaml`** – Main configuration file specifying pipeline parameters, input/output locations, and reference genome details.
-   **`fastq_paths.txt`** – Lists sample IDs and their corresponding FASTQ file paths for processing.

## Configuration Details

The `config.yaml` file includes settings for:

**General Settings:**

-   **`threads`**: The number of threads to use for parallel execution of pipeline steps (default: 36).
-   **`output_path`**: The main directory where pipeline outputs will be stored (default: `/scratch/indikar_root/indikar1/shared_data/ont_10x_pipeline_test/`).

## Input Files

- **`inputs`**:
  - **`fastq_file_paths`**: Path to the `fastq_paths.txt` file, which contains a list of sample IDs and their corresponding FASTQ file paths. This file is crucial for linking sample identifiers to their sequencing data.
  - **`inputs_are_directories`**:  **(Located in `config.yaml`)** A boolean value (`true` or `false`) that specifies whether the paths in `fastq_paths.txt` refer to individual FASTQ files or directories containing FASTQ files. 
  - **`reference_fasta`**: Path to the reference genome FASTA file used for alignment (default: `/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.cdna.all.fa`).
  - **`reference_annotation`**: Path to the reference genome annotation GTF file used for gene quantification and other downstream analyses (default: `/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.107.gtf`).

### FASTQ File Merging

- **Pros:**
  1. Merging allows the user to control the size and number of files.
  2. Facilitates combining data from multiple sequencing runs.

- **Cons:**
  1. Automated merging assumes all files in the directory should be merged.
  2. The directory should contain **ONLY** FASTQ files to avoid errors.


**Demultiplexing Parameters (using `demux`):**

-   **`expected_cells`**: The anticipated number of cells in the sample (default: 10000). This helps optimize the demultiplexing process.
-   **`barcode_whitelist`**: Path to a file containing the list of valid cell barcodes (default: `/nfs/turbo/umms-indikar/shared/projects/cell_cycle/data/10x_barcode_lists/3M-3pgex-may-2023.txt`). This file is used to filter out reads with incorrect or unknown barcodes.

**Alignment Parameters (using `minimap2`):**

-   **`alignment`**:
    -   **`minimap2_params`**: Parameters passed directly to the `minimap2` aligner, controlling alignment behavior. The default setting `"-ax splice -uf --secondary=no --MD"` indicates settings optimized for spliced alignment of RNA-seq data, including consideration for strand information and suppression of secondary alignments (default: `-ax splice -uf --secondary=no --MD`).
    -   **`mapq_threshold`**: Minimum mapping quality (MAPQ) score for a read to be considered aligned (default: 1). Reads below this threshold are discarded.

**Gene Quantification Parameters**

-   **`counts`**:
    -   **`umi_distance`**: The maximum edit distance allowed between UMIs to be considered the same during gene quantification. (default: 1).

**`fastq_paths.txt` File Format:**

This file should be a tab-separated or space-separated text file with two columns:

1.  **Sample ID:** A unique identifier for each sample.
2.  **FASTQ File Path:** The absolute path to the corresponding FASTQ file for that sample.

## SLURM Configuration

To run the pipeline on a SLURM-managed cluster, you need the following configuration files:

### `cluster.json`
Defines job-specific parameters for each sub-job in the workflow, including:
- **Logging directory** – Path for storing job logs  
- **SLURM user account** – The user account under which jobs will run  
- **SLURM billing account** – The account used for job billing  

### `cluster/config.yaml`
Specifies the default SLURM submission parameters for sub-jobs. This file typically does not require modification.

Ensure both files are correctly configured before executing the workflow.
