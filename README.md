# ONT 10x Feature Barcoding Pipeline

This pipeline is designed for processing single-cell transcriptomics data generated using Oxford Nanopore Technologies (ONT) sequencing with 10x Genomics GEM chips. It focuses on disambiguating BioLegend feature barcodes.

## Overview

The pipeline is tailored to handle the unique characteristics of ONT single-cell transcriptomics data prepared using 10x Genomics technology. It ensures efficient processing, barcode disambiguation, and downstream analysis readiness.

## Cloning the Pipeline Repository

To obtain the pipeline, clone the repository from GitHub:

```bash
git clone https://github.com/CooperStansbury/ont_10x_feature_barcodes.git
cd ont_10x_feature_barcodes
```

This ensures you have the latest version of the pipeline.

## Running the Pipeline

To run the pipeline, follow these steps:

1.  **Set Up the Environment**  
    Ensure Conda is installed and set up correctly. The pipeline requires the following Conda environment:

    *   **Top-level environment:** `workflow-env.yaml`

    Install the environments by running:

    ```bash
    conda env create -f envs/workflow-env.yaml
    ```

    After installation, activate the top-level environment:

    ```bash
    conda activate workflow-env
    ```

2.  **Prepare Configuration Files**  
    Update the `config.yaml` file with the correct paths to your input data, reference genome, and output directories.

3.  **Execute the Pipeline**  
    Always run Snakemake using the `--use-conda` flag to ensure proper dependency management:

    ```bash
    snakemake --use-conda --cores <num_cores>
    ```

    Replace `<num_cores>` with the number of threads available for computation.

4.  **Dry Run (Optional)**  
    To verify the workflow without executing commands:

    ```bash
    snakemake --use-conda --configfile config.yaml -n
    ```

5.  **Cluster Execution (Optional)**  
    If running on an HPC system, submit jobs using:

    ```bash
    snakemake --use-conda --cluster "sbatch --mem={resources.mem_mb}" --jobs 10
    ```

## Input Requirements
TODO

## Output
TODO

## Troubleshooting

*   **Ensure Conda is installed and environments are set up correctly.**
*   **Verify input file paths:** Ensure all paths in `config.yaml` and `fastq_paths.txt` are correct.
*   **Check for missing dependencies:** Run `snakemake --use-conda --conda-create-envs-only`.
*   **Review logs:** Check the output logs in the `logs` directory for error messages and troubleshooting hints.

---

For further details and support, refer to the official Snakemake documentation:  
[https://snakemake.readthedocs.io](https://snakemake.readthedocs.io)