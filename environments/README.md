# Conda Environments

This directory contains YAML configuration files for setting up various Conda environments. Each file can be used to create a Conda environment with the specified dependencies and configurations.

## Usage

### Creating a New Environment

To create a Conda environment from a YAML file, use the following command:

```bash
conda env create -f environment.yml
```

Replace `environment.yml` with the appropriate YAML file for your desired environment.

### Updating an Existing Environment

To update an existing Conda environment with changes in the YAML file, run:

```bash
conda env update --file environment.yml --prune
```

The `--prune` option removes packages that are no longer required.

### Activating an Environment

Once the environment is created, activate it using:

```bash
conda activate my_env
```

Replace `my_env` with the environment name specified in the YAML file.

### Exporting an Environment

To export an existing Conda environment to a YAML file for sharing or backup purposes:

```bash
conda env export > environment.yml
```

If you wish to exclude specific system-specific dependencies (like paths), add:

```bash
conda env export --no-builds > environment.yml
```

### Removing an Environment

To remove an environment that is no longer needed:

```bash
conda env remove -n my_env
```

### Listing Available Environments

To view all Conda environments on your system:

```bash
conda env list
```

Or:

```bash
conda info --envs
```

## File Structure

This folder contains the following environment YAML files:

- `environment.yml` – General environment with core dependencies.
- `dev_environment.yml` – Development environment with additional tools for debugging and testing.
- `ml_environment.yml` – Environment tailored for machine learning projects.
- `data_analysis.yml` – Environment for data analysis tasks with relevant libraries.

## Notes

- Ensure Conda is installed and accessible from your terminal.
- Some environments may require specific versions of Conda or Python.
- Consider using `mamba` instead of `conda` for faster dependency resolution by replacing `conda` with `mamba` in the above commands.

## Troubleshooting

If you encounter issues during environment installation or activation, consider the following solutions:

1. **Resolve Dependency Conflicts:**  
   Run the following command to check for dependency conflicts:  
   ```bash
   conda config --set channel_priority flexible
   ```

2. **Update Conda:**  
   Ensure Conda is up-to-date:  
   ```bash
   conda update conda
   ```

3. **Remove and Recreate the Environment:**  
   If issues persist, remove the environment and recreate it:  
   ```bash
   conda env remove -n my_env
   conda env create -f environment.yml
   ```

4. **Check Conda Configuration:**  
   View your Conda configuration settings using:  
   ```bash
   conda config --show
   ```

## Additional Resources

- Official Conda documentation: [https://docs.conda.io](https://docs.conda.io)
- Mamba documentation: [https://mamba.readthedocs.io](https://mamba.readthedocs.io)

---

Feel free to modify or add more environment files to suit your needs.
