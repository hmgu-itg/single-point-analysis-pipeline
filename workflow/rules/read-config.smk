import sys
import pandas as pd

configfile: "config.yaml"


if isinstance(config['phenotypes'], str):
    config['phenotypes'] = pd.read_csv(config['phenotypes'], header = None)[0].to_list()
