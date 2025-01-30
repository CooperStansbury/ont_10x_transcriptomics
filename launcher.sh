#!/bin/bash

#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --time=36:00:00

# Default values
CONFIG='config/cluster'
CORES=36

# Parse keyword arguments with flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--config)
      CONFIG="$2"
      shift # past argument
      shift # past value
      ;;
    -j|--cores)
      CORES="$2"
      shift # past argument
      shift # past value
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Print argument values nicely
echo "---------------------------------"
echo "  Snakemake Workflow Settings"
echo "---------------------------------"
echo "Config File: ${CONFIG}"
echo "Number of Cores: ${CORES}"
echo "---------------------------------"

## build the workflow from the most current snakefile
cp Snakefile workflow.smk

# run it (Uncomment to run)
snakemake --profile ${CONFIG} --use-conda --cores ${CORES} --rerun-incomplete --latency-wait 90 --verbose -s workflow.smk