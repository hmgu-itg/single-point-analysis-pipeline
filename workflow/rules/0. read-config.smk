"""
Snakefile 0.


"""
import sys
from pathlib import Path

configfile: "config.yaml"

BFILE = config['bfile']
BFILE_INPUTS = multiext(BFILE, '.bed', '.bim', '.fam')
GRM = config['grm']
PHENOTYPE_FILE = config['phenotype_file']

OUTPUT_PATH = config['output_path']

LMISS = config['lmiss']

config['p-value'] = float(config['p-value'])
