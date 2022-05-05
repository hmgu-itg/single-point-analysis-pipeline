"""
Snakefile 0.


"""
import sys
import pandas as pd

configfile: "config.yaml"


if isinstance(config['phenotypes'], str):
    config['phenotypes'] = pd.read_csv(config['phenotypes'], header = None)[0].to_list()


config['QC_thresholds']['MAC'] = int(config['QC_thresholds']['MAC'])
config['QC_thresholds']['HWE'] = float(config['QC_thresholds']['HWE'])
config['QC_thresholds']['missingness'] = float(config['QC_thresholds']['missingness'])
config['QC_thresholds']['p-value'] = float(config['QC_thresholds']['p-value'])
