import pandas as pd
import numpy as np
import os
import sys
import pysam

def process_bam(in_bam, out_bam, out_record):
    """
    Processes a BAM file to extract barcode, UMI, and read name information, 
    and stores this information in a separate CSV file.

    This function iterates through each alignment in the input BAM file, extracts
    the barcode, UMI, and read name from the query name, and adds these as tags
    to the alignment record. The modified alignment is then written to a new BAM file.
    Additionally, relevant information about each alignment is stored in a list of
    dictionaries, which is ultimately converted to a pandas DataFrame and saved as a CSV file.

    Args:
        in_bam (str): Path to the input BAM file.
        out_bam (str): Path to the output BAM file.
        out_record (str): Path to the output CSV file.
    """
    records = []

    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in:
                barcode = align.query_name.split('_')[0]
                umi = align.query_name.split('_')[1].split("#")[0]
                read_name = align.query_name.split("#")[1][:-2]

                # Add tags to the alignment record
                align.set_tag('CB', barcode, value_type="Z")
                align.set_tag('UB', umi, value_type="Z")
                align.set_tag('RD', read_name, value_type="Z")
                bam_out.write(align)

                # Structure the record table
                row = {
                    'barcode' : barcode,
                    'umi' : umi,
                    'read_name' : read_name,
                    'forward' : align.is_forward,
                    'mapping_quality' : align.mapping_quality,
                    'query_length' : align.query_length,
                }

                records.append(row)

    # Save the record table
    records = pd.DataFrame(records)
    records.to_csv(out_record)

if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    out_record = sys.argv[3]

    process_bam(in_bam, out_bam, out_record)