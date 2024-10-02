#!/usr/bin/env python3

__version__ = '0.2'

import sys

import pandas as pd
from peakplotter.peakit import peakit
from peakplotter.plotpeaks import read_assoc, get_signals
from peakplotter.test_utils import get_test_logger

assocfile = sys.argv[1]
flank_bp = int(sys.argv[2])
signif = float(sys.argv[3])
output = sys.argv[4]


chr_col = 'Chr'
pos_col = 'bp'
rs_col = 'SNP'
pval_col = 'p'
a1_col = 'A1'
a2_col = 'A2'
maf_col = 'Freq' # Frequency of A1
logger = get_test_logger()


print(f'Running collect_peaks.py (v{__version__})')
print(f'''
chr_col = '{chr_col}'
pos_col = '{pos_col}'
rs_col = '{rs_col}'
pval_col = '{pval_col}'
a1_col = '{a1_col}'
a2_col = '{a2_col}'
maf_col = '{maf_col}'
''')

assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
print('Getting signals')
signals = get_signals(assoc, signif, chr_col, pos_col, pval_col)
print(signals)
if signals.empty:
    print("No peaks found. Exiting.")
    pd.DataFrame().to_csv(output, sep = '\t', header = False, index = False)
    sys.exit(0)
print('Making peak_collections')
peak_collections = peakit(signals, pval_col, chr_col, pos_col, flank_bp)
print('peak_collections before merge')
print(peak_collections)
peak_collections.merge()
print('peak_collections after merge')
print(peak_collections)
peaked = peak_collections.data
peaked.to_csv(output, sep = '\t', header = False, index = False)