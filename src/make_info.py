# This script create .info.gz file from bgzipped vcf file.
# info files are already included if all output files are downloaded from TOPmed imputation server
# Requires bcftools
import argparse
import os
from xopen import PipedCompressionWriter, xopen
lst_arg = ['--input', '--output_dir', '--output_fn', '--verbose'] # Available arguments
parser = argparse.ArgumentParser(description='Generate info files')
parser.add_argument('--input', nargs='*', type=str, required=True,
                    help='(Required) Multiple input files are allowed. Must in gzipped or bgziped VCF format')
parser.add_argument('--output_dir', type=str, default='',
                    help='(Optional) Directory for output files. Default is current working directory.')
parser.add_argument('--output_fn', type=str, nargs='*',
                    help='(Optional) Default is input file name with suffix replaced by ".info.gz".')
parser.add_argument('--col_names', nargs='*', default=['AF', 'MAF', 'R2', 'IMPUTED/TYPED/TYPED_ONLY'],
                    help="(Optional) Default is ['AF', 'MAF', 'R2', 'IMPUTED/TYPED/TYPED_ONLY']. Column names of alt frequency, MAF, imputation quality score, genotyped. Separated by space")
parser.add_argument('--thread', default=1, type=int,
                    help='(Optional) Default is 1. Defines how many thread to use in multiprocessing. If number of threads <0, will use 1 instead of user supplied value.')
parser.add_argument('--write_with', default='bgzip',
                    help='(Optional) User specified program to write compressed file. Default is bgzip. (gzip can be used unless required to bgzip info files)')
parser.add_argument('--verbose', default=False,
                    help='(Optional) Print help message if True. Default is False. Valid values are (not case-sensitive): {0|1|True|False}')
args = parser.parse_args()
if (not args.verbose) or args.verbose.upper()=='FALSE' or args.verbose=='0':
    args.verbose = False
else:
    args.verbose = True
    parser.print_help()  # Print help message if desired

if args.output_fn is None: # Fill output file names as input name + '.info.gz'
    args.output_fn = [os.path.split(fn)[1].split('.')[0] + '.info.gz' for fn in args.input]
if args.thread<0: args.thread = 1

print('\nOptions used')
for arg in lst_arg:
    print('\t' + arg, eval('args.'+arg[2:]))
print('\n')

# ################################ Sanity checks ################################
if args.verbose:
    print('#### Sanity check: checking user inputs')

invalid_flag = False # Track if any input value is invalid
# Check if input files exist
for fn in args.input:
    if not os.path.isfile(fn):
        print(f'Error: Invalid input file name: {fn} is not a file or does not exist')
        invalid_flag = True

# Check if number of input and output files matches
if len(args.input) != len(args.output_fn):
    print(f'Error: Number of input file ({len(args.input)}) does not match number of output file ({len(args.output_fn)})')
    invalid_flag = True

if len(args.col_names) > 4:
    print(f'Error: Too many values. Only need (4) column names of alt frequency, MAF, imputation quality score and genotyped.')
    invalid_flag = True
if len(args.col_names) < 4:
    print(f'Error: Not enough values. Need (4) column names of alt frequency, MAF, imputation quality score and genotyped. Values should be separated by space.')
    invalid_flag = True

if invalid_flag:
    print('Exit')
    exit()

if args.verbose:
    print('#### Sanity check: Done. All inputs are valid')
# ################################ End of sanity check ################################

# ################################ Write to output file ################################
# Create directory if user provide --output_dir does not exist
if len(args.output_dir) and (not os.path.isdir(args.output_dir)):
    if args.verbose:
        print(f'#### Output directory {args.output_dir} does not exist. Creating output directory')
    os.mkdir(args.output_dir)

for i in range(len(args.input)):
    if args.verbose:
        print(f'\n#### Processing files: Start. input file = {args.input[i]}; output file = {os.path.join(args.output_dir, args.output_fn[i])}')

    output_fh = PipedCompressionWriter(path=os.path.join(args.output_dir + args.output_fn[i]), threads_flag="-@", program_args=[args.write_with],
                           threads=args.thread)
    # Need these columns: 'SNP', 'REF(0)', 'ALT(1)', 'Genotyped','ALT_Frq', 'MAF', 'Rsq'
    output_fh.write('SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tRsq\tGenotyped\n')
    if args.verbose:
        print(f'#### Processing files: Write column headers to output file {args.output_fn[i]}')

    with xopen(args.input[i]) as input_fh:
        line = input_fh.readline().strip()
        while line[:2] == '##': # Read through header lines start with ##
            line = input_fh.readline().strip()
        column_names = line.split()
        # Find index of desired columns
        idx_snp = column_names.index('ID')
        idx_ref = column_names.index('REF')
        idx_alt = column_names.index('ALT')
        idx_info = column_names.index('INFO') # INFO column

        # Get indices of AF (alt_frq), MAF, Rsq and genotyped from INFO column
        line = input_fh.readline().strip()  # Read lines that contain actual values
        info_val = line.split()[idx_info]
        idx = 0 # Track index in INFO column (fields separated by ';')
        idx_alt_frq, idx_maf, idx_rsq, idx_genotyped = -1, -1, -1, -1 # FInd indices of these 4 columns
        for val in info_val.split(';'):
            key = val.split('=')[0]
            if key == args.col_names[0]: # AF
                idx_alt_frq = idx
            elif key == args.col_names[1]: # MAF
                idx_maf = idx
            elif key == args.col_names[2]: # R2
                idx_rsq = idx
            elif key in args.col_names[3].split('/'): # IMPUTED/TYPED/TYPED_ONLY
                idx_genotyped = idx
            idx += 1
        # Sanity check: exit if any field is missing
        if idx_alt_frq==-1 or idx_maf==-1 or idx_rsq==-1 or idx_genotyped==-1:
            print(f"Error: Missing field(s) indicated by -1 (AF, MAF, Rsq, Genotyped): {idx_alt_frq}, {idx_maf}, {idx_rsq}, {idx_genotyped}\nExit")
            exit()

        if args.verbose:
            print(f'#### Processing files: Indices of columns (name (index)): ID ({idx_snp}); REF ({idx_ref}); ALT ({idx_alt}); INFO ({idx_info})')
            print(f'#### Processing files: Indices of fields in INFO column (name (default name) (index)): {args.col_names[0]} (AF) ({idx_alt_frq}); {args.col_names[1]} (MAF) ({idx_maf}); {args.col_names[2]} (R2) ({idx_rsq}); {args.col_names[3]} (IMPUTED/TYPED/TYPED_ONLY) ({idx_genotyped})')
            print(f'#### Processing files: Loop through input file')

        count = 0
        while line != '':
            tmp_lst = line.split()
            snp = tmp_lst[idx_snp]
            ref = tmp_lst[idx_ref]
            alt = tmp_lst[idx_alt]
            info_val = tmp_lst[idx_info].split(';')
            alt_frq = info_val[idx_alt_frq].split('=')[1]
            maf = info_val[idx_maf].split('=')[1]
            rsq = info_val[idx_rsq].split('=')[1]
            genotyped = info_val[idx_genotyped]
            output_fh.write(f'{snp}\t{ref}\t{alt}\t{alt_frq}\t{maf}\t{rsq}\t{genotyped}\n')
            count += 1
            if args.verbose:
                if count%1000 == 0:
                    print('.', end='', flush=True)
                if count%5000 == 0:
                    print(f'{count} lines processed')
            line = input_fh.readline().strip()
    output_fh.close()
    if args.verbose:
        print(f'#### Processing files: Done. input={args.input[i]}')

