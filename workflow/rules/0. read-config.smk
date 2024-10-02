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

LMISS = config['lmiss']

config['p-value'] = float(config['p-value'])
