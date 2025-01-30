# Executing on UM Great Lakes

This guide provides the necessary steps to configure and execute jobs on the University of Michigan's Great Lakes HPC cluster.

## Steps

1. **Update the Cluster Configuration**  
   Modify `config/cluster.json` to specify per-job parameters, including SLURM account and user information.

2. **Update the Job Launcher Script**  
   Ensure the SLURM account and user details are correctly set in the header of `launcher.sh`.

3. **Check `.bashrc`**  
   Verify that your environment is properly set up for execution, including necessary module loads and paths.

## Additional Resources

- [SchlossLab: Great Lakes SLURM Guide](https://github.com/SchlossLab/Great_Lakes_SLURM)  
- [Kelly Sovacool: Snakemake HPC Minimum Working Example](https://github.com/kelly-sovacool/snakemake_hpc_mwe)
