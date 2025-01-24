# Scripts

This directory contains utility scripts for processing genomic data.

## `compile_anndata.py`

Combines multiple `h5ad` files (e.g., from different chromosomes), merges gene metadata from a table, and saves the result as a single `h5ad` file.

## `tag_bam.py`

Processes a BAM file to extract barcode, UMI, and read name information from the query name. Adds this information as tags to the alignment record and saves it to a new BAM file. Also, it generates a CSV file containing relevant alignment details.