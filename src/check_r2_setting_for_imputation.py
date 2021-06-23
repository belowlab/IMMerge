# This code is to check what r2 setting was applied for TOPMed imputation
# This may not be necessary if merging allows missing values
# Only need to run once

import gzip
import sys
from process_args import process_args

# File name can be passed to this code in terminal, or use import this code as in a script (need to modify a little)
args = sys.argv
verbose=True

# Process terminal input
# dict_flags contains values for below flags:
#   --input, --output, --verbose, --missing, --r2_threshold, --r2_output
dict_flags = process_args(args) # Process terminal input

if dict_flags['--verbose']!='true': verbose=False

# ----------------------- Helper functions -----------------------
# This function check imputation setting info based on header lines of input .dose.vcf.gz files
# Parameter: fn: input file name
# Output: print imputation settings of each input file to console
def check_imputatilson_parameters(lst_fn=dict_flags['--input']):
    print('\nCheck r2 filter used for imputation:')
    try:
        # Read in file headers of each input file
        for fn in lst_fn:
            fh = gzip.open(fn, 'rt')
            line = fh.readline()
            # Header lines start with '##'
            while line[0:2]=='##':
                if line[0:4]=='##r2':
                    print('\t', fn.split('/')[-1]+':', line.strip()[2:])
                line = fh.readline()
            fh.close()
    except:
        print('Error: File not found:', fn)
        raise IOError('Input file not found')

check_imputatilson_parameters()
