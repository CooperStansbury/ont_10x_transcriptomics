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


HEADER_STR = "#" * 20


def log(message):
    """Prints a timestamped message to the console."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp} - {message}")


def check_consistent_extensions(basenames):
  """
  Checks if a list of basenames have consistent extensions.

  Args:
    basenames: A list of file basenames (e.g., ["file1.txt", "file2.txt.gz"]).

  Returns:
    The consistent extension (including any leading periods) if found.

  Raises:
    ValueError: If the extensions are not consistent or if no files are provided.
  """
  if not basenames:
    raise ValueError("No basenames provided.")

  extensions = []
  for basename in basenames:
    parts = basename.split(".")
    if len(parts) > 2 and parts[-2] == "fastq" and parts[-1] == "gz":
      extensions.append(".fastq.gz")
    else:
      extensions.append(os.path.splitext(basename)[1])

  extensions_set = set(extensions)
  if len(extensions_set) > 1:
    raise ValueError(f"Inconsistent extensions found: {', '.join(extensions_set)}")

  return extensions_set.pop()


def get_output_filenames(input_df, output_path):
  """
  Generates output filenames by combining sample IDs from an input DataFrame 
  with a consistent file extension and the output directory path.

  Args:
    input_df: A pandas DataFrame with 'sample_id' and 'file_path' columns.
    output_path: The directory where output files will be written.

  Returns:
    A list of complete output filenames (including the output path and extension).

  Raises:
    ValueError: If the 'file_path' column contains inconsistent extensions.
  """
  output_names = input_df['sample_id'].to_list()
  basenames = input_df['file_path'].apply(os.path.basename).to_list()
  extension = check_consistent_extensions(basenames)

  output_names = [f"{output_path}fastq/{x}{extension}" for x in output_names]
  return output_names


def read_names(names_file):
    """
    Reads a file containing names, one per line.

    Args:
        names_file (str): Path to the file containing names.

    Returns:
        list: A list of names.
    """
    with open(names_file) as f:
        return [line.strip() for line in f]
    
    
    
    
    