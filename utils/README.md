# pipeline_utils.py

This file contains helper functions for your data pipelines.

## Functions

* **log(message):** Prints a message with a timestamp.
* **check_consistent_extensions(basenames):** Checks if a list of file basenames have consistent extensions.
* **get_output_filenames(input_df, output_path):** Generates output filenames by combining sample IDs from a DataFrame and a consistent file extension.