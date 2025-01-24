import pandas as pd
import numpy as np
import os



if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    out_record = sys.argv[3]

    process_bam(in_bam, out_bam, out_record)