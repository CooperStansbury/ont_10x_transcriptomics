# ONT 10x Transcriptomics Pipeline

This pipeline processes single-cell transcriptomics data generated using Oxford Nanopore Technologies (ONT) sequencing with 10x Genomics GEM chips.

## Overview

The pipeline assumes that the sequencing library was prepared using 10x Genomics technology and is designed to handle the specific requirements of ONT single-cell transcriptomics data.

## Running the Pipeline

To run the pipeline, follow these steps:

1. **Set Up the Environment**  
   Ensure Conda is installed and set up correctly. The pipeline requires the following Conda environments:
   
   - **Top-level environment:** `workflow-env.yaml`
   - **Core pipeline environment:** `pipeline-core.yaml`
   - Additional specific environments for different pipeline components.

   Install and activate the environments by running:

   ```bash
   conda env create -f workflow-env.yaml
   conda env create -f envs/pipeline-core.yaml
   ```

2. **Prepare Configuration Files**  
   Update the `config.yaml` file with the correct paths to your input data, reference genome, and output directories.

3. **Execute the Pipeline**  
   Always run Snakemake using the `--use-conda` flag to ensure proper dependency management:

   ```bash
   snakemake --use-conda --configfile config.yaml --cores <num_cores>
   ```

   Replace `<num_cores>` with the number of threads available for computation.

4. **Dry Run (Optional)**  
   To verify the workflow without executing commands:

   ```bash
   snakemake --use-conda --configfile config.yaml -n
   ```

5. **Cluster Execution (Optional)**  
   If running on an HPC system, submit jobs using:

   ```bash
   snakemake --use-conda --cluster "sbatch --mem={resources.mem_mb}" --jobs 10
   ```

## Input Requirements

- **Raw FASTQ Files:** Path specified in `config.yaml`
- **Reference Genome:** FASTA and GTF files for alignment and annotation
- **Configuration File:** Contains parameters for alignment, filtering, and output directories

## Output

The pipeline generates the following outputs:

- **Aligned Reads:** Processed and aligned sequences
- **Gene Expression Matrices:** Quantified expression levels
- **Quality Control Reports:** Metrics and visualizations

## Troubleshooting

- **Ensure Conda is installed and environments are set up correctly.**
- **Verify input file paths:** Ensure all paths in `config.yaml` are correct.
- **Check for missing dependencies:** Run `snakemake --use-conda --conda-create-envs-only`.
- **Review logs:** Check the output logs for error messages and troubleshooting hints.

---

For further details and support, refer to the official Snakemake documentation:  
[https://snakemake.readthedocs.io](https://snakemake.readthedocs.io)
