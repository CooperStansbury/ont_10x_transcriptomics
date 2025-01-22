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
INPUT_ABS_PATH = os.path.abspath(config['inputs'])
INPUT_BASENAME = os.path.basename(INPUT_ABS_PATH)
input_df = pd.read_csv(INPUT_ABS_PATH, comment="#")
samples = input_df['sample_id'].to_list()

print(f"\n{pu.HEADER_STR} FEATURE VALUES {pu.HEADER_STR}")

print(
    tabulate.tabulate(
        input_df, 
        headers='keys', 
        tablefmt='plain',
        showindex=False,
    )
)


rule all:
    input:
        OUTPUT_PATH + "config/" + CONFIG_BASENAME,
        OUTPUT_PATH + "config/" + INPUT_BASENAME,


rule get_config:
    input:
        CONFIG_ABS_PATH
    output:
        OUTPUT_PATH + "config/" + CONFIG_BASENAME,
    shell:
        """ cp {input} {output} """


rule get_input_list:
    input:
        INPUT_ABS_PATH
    output:
        OUTPUT_PATH + "config/" + INPUT_BASENAME,
    shell:
        """ cp {input} {output} """
    
        


