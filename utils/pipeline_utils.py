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