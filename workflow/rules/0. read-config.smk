"""
Snakefile 0.


"""
import sys
from pathlib import Path

configfile: "config.yaml"

COHORT = config['cohort']
GROUP = config['group']
BFILE = config['bfile']
GRM = config['grm']
PHENOTYPE = config['phenotype']
PHENOTYPE_FILE = config['phenotype_file']

# Checks
if not isinstance(COHORT, str) or len(COHORT)<2:
    print('`cohort` config value must be string and length>1')
    sys.exit(1)

if not isinstance(config['group'], str) or len(COHORT)<2:
    print('`group` config value must be string and length>1')
    sys.exit(1)


config['QC_thresholds']['MAC'] = int(config['QC_thresholds']['MAC'])
config['QC_thresholds']['HWE'] = float(config['QC_thresholds']['HWE'])
config['QC_thresholds']['missingness'] = float(config['QC_thresholds']['missingness'])
config['QC_thresholds']['p-value'] = float(config['QC_thresholds']['p-value'])
