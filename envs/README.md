# Conda Environments

This directory contains YAML configuration files for setting up Conda environments required to run the pipeline.

## Environment Overview

- **`workflow-env.yaml`** – Required for executing the Snakemake pipeline.
- **`pipeline-core.yaml`** – Contains software dependencies needed to run individual Snakemake rules.

## Setup Instructions

Before running the pipeline, create the necessary Conda environments:

```bash
conda env create -f <env.yaml>
```

Replace `<env.yaml>` with the specific environment file (e.g., `workflow-env.yaml`).

## Common Conda Commands

### Updating an Environment

To update an existing environment based on changes in the YAML file:

```bash
conda env update --file <env.yaml> --prune
```

The `--prune` option removes unnecessary packages.

### Activating an Environment

After creation, activate the environment with:

```bash
conda activate <env_name>
```

Replace `<env_name>` with the environment name specified in the YAML file.

### Listing Available Environments

To check available Conda environments:

```bash
conda env list
```

## Additional Resources

- [Conda Documentation](https://docs.conda.io)
- [Mamba Documentation](https://mamba.readthedocs.io)

---

Modify or add environment files as needed to fit your project requirements.
