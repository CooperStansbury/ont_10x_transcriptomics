import os
import sys
import glob
import re
from datetime import datetime
from pathlib import Path
import pandas as pd
import yaml
import json
import tabulate 

""" PATH CONFIG """
BASE_DIRECTORY = Path(workflow.basedir)

# config details
CONFIG_PATH = "/config/config.yaml"
CONFIG_BASENAME = os.path.basename(CONFIG_PATH)
CONFIG_ABS_PATH = str(BASE_DIRECTORY) + CONFIG_PATH
configfile: CONFIG_ABS_PATH 

# pipeline utilities
UTILS_PATH = "/utils/"
sys.path.append(str(BASE_DIRECTORY) + UTILS_PATH) 
import pipeline_utils as pu

""" PRINT THE EXECUTION DETAILS """
print(f"\n{pu.HEADER_STR} EXECUTION DETAILS {pu.HEADER_STR}")
pu.log(f"Pipeline started")
print(f"Base directory: {BASE_DIRECTORY}")
print(f"Config file path: {CONFIG_ABS_PATH}")

# Print config values
print(f"\n{pu.HEADER_STR} CONFIG DETAILS {pu.HEADER_STR}")
print(json.dumps(config, indent=4)) 

""" HELPER VARIABLES """
OUTPUT_PATH = config['output_path']

""" LOAD INPUTS """
INPUT_ABS_PATH = os.path.abspath(config['inputs']['fastq_file_paths'])
INPUT_BASENAME = os.path.basename(INPUT_ABS_PATH)
input_df = pd.read_csv(INPUT_ABS_PATH, comment="#")
samples = input_df['sample_id'].to_list()

# get new path names
input_file_paths = input_df['file_path'].to_list()
output_file_paths = pu.get_output_filenames(input_df, OUTPUT_PATH)
extension = pu.check_consistent_extensions(input_file_paths)

print(f"\n{pu.HEADER_STR} INPUT FILES {pu.HEADER_STR}")
for _, row in input_df.iterrows():
    fbasename = os.path.basename(row['file_path'])
    print(f"{row['sample_id']}: {fbasename} ({row['file_path']})")


""" RULE FILES """
include: "rules/gather.smk"
include: "rules/pipeline-core.smk"


rule all:
    input:
        OUTPUT_PATH + "config/" + CONFIG_BASENAME,
        OUTPUT_PATH + "config/" + INPUT_BASENAME,
        OUTPUT_PATH + 'references/genome_build.txt',
        OUTPUT_PATH + 'references/reference.mmi',
        OUTPUT_PATH + "references/barcode_whitelist.txt",
        output_file_paths,
        expand(OUTPUT_PATH + "demultiplex/{sid}.done", sid=samples),
        OUTPUT_PATH + "reports/seqkit_stats/raw_fastq_report.txt",
        OUTPUT_PATH + 'reports/seqkit_stats/demultiplexed_fastq_report.txt',
        expand(OUTPUT_PATH + "mapping/{sid}.tagged.bam.bai", sid=samples),



