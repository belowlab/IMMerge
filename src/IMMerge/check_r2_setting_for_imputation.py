# Step 00 (optional, set to default in this version)
# This code is to check what r2 setting was applied for TOPMed imputation
# This may not be necessary if merging allows missing values
# Only need to run once
import gzip

# This function check imputation setting info based on header lines of input .dose.vcf.gz files
# Parameter: fn: input file name (stored in dict_flags['--input'])
# Output:
# - Print imputation settings of each input file to console
# - Processed arguments are returned in dict_flags
def check_imputatilson_parameters(lst_fn):
    print('\nCheck r2 filter used for imputation:')
    try:
        # Read in file headers of each input file
        for fn in lst_fn:
            fh = gzip.open(fn, 'rt')
            line = fh.readline()
            # Header lines start with '##'
            while line[0:2]=='##':
                if line[0:4]=='##r2':
                    print('\t'+fn.split('/')[-1]+':', line.strip()[2:])
                line = fh.readline()
            fh.close()
        # return dict_flags
    except:
        print('Error: File not found:', fn, '\n')
        raise IOError('Input file not found')

# check_imputatilson_parameters()
