# Conda Environments

This directory contains YAML configuration files for setting up Conda environments required to run the pipeline.

## Environment Overview

-   **`workflow-env.yaml`**  –  **Required:** Contains dependencies for executing the Snakemake workflow manager. This is the only environment you need to create manually.

**Automatically Managed Environments (Created by Snakemake):**

-   **`pipeline-core.yaml`**  – Contains core software dependencies needed for individual Snakemake rules.
-   **`pipeline-QC.yaml`**  – Contains software dependencies needed for Snakemake quality control and reporting.

**Important:** You do  **not**  need to create or activate the `pipeline-core.yaml` or `pipeline-QC.yaml` environments manually. Snakemake automatically creates and manages these environments based on their respective YAML files when you run the pipeline.

## Setup Instructions

Before running the pipeline, you  **only need to create the `workflow-env`**:


```bash
conda env create -f workflow-env.yaml
```

## Common Conda Commands

### Activating an Environment

After creation, activate the environment with:

```bash
conda activate <env_name>
```

Replace `<env_name>` with the environment name specified in the YAML file.


## Additional Resources

- [Conda Documentation](https://docs.conda.io)
- [Mamba Documentation](https://mamba.readthedocs.io)

---

Modify or add environment files as needed to fit your project requirements.
