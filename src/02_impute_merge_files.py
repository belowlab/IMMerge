# This code is used to merge imputed files from TopMed based on list of saved variants in prior steps
# The goal is to combine post-imputation vcf.gz files into one vcf.gz file



print('\nVersion 1.0')

import gzip
import sys
from process_args import process_args
import get_SNP_list
import check_r2_setting_for_imputation

# File name can be passed to this code in terminal, or use import this code as in a script (need to modify a little)
args = sys.argv
verbose=True
# Process terminal input
# dict_flags contains values for below flags:
#   --input, --output, --verbose, --missing, --r2_threshold, --r2_output
dict_flags = process_args(args) # Process terminal input
if dict_flags['--verbose']!='true': verbose=False

check_r2_setting_for_imputation.check_imputatilson_parameters(lst_fn=dict_flags['--input'])
get_SNP_list.get_snp_list(dict_flags)


# ----------------------- Helper functions -----------------------

# This function reads in a file with one snp in each line (ie. variants_kept.txt),
# then return a list of all SNPs.
# The first line needs to be skipped (header line)
def get_snp_lst(variants_kept_fn='variants_kept.txt'):
    lst_snp = []
    with open(variants_kept_fn, 'r') as fh:
        line = fh.readline()  # Skip the header line
        line = fh.readline().strip()
        while line != '':
            lst_snp.append(line.split()[0])
            line = fh.readline().strip()
    return lst_snp

# Merge input files together
def merge_files():
    print('\nStart merging files of:')
    lst_fn = dict_flags['--input']
    output_fn = dict_flags['--output'] + '.dose.vcf.gz'

    # Read in SNPs need to be kept from file
    lst_snp_to_keep = get_snp_lst()

    fh_group1 = gzip.open(input_dir + lst_fn[0], 'rt')
    fh_group2 = gzip.open(input_dir + lst_fn[1], 'rt')
    fh_group3 = gzip.open(input_dir + lst_fn[2], 'rt')
    fh_output = gzip.open(output_fn, 'wt')  # Open for writing (appending)

    # Skip the first 18 lines (header lines of imputation info).
    # Real data starts from line 19, but need to write the header lines into merged file once
    for i in range(18):
        line_group1 = fh_group1.readline().strip()
        line_group2 = fh_group2.readline().strip()
        line_group3 = fh_group3.readline().strip()
        fh_output.write(line_group3 + '\n')

    # Need to write the header line of columns into file
    # Need to merge group1, group2 and group3 column headers
    line_group1 = fh_group1.readline().strip()
    line_group2 = fh_group2.readline().strip()
    line_group3 = fh_group3.readline().strip()
    line_to_write = line_group1 + '\t' + line_group2.split(maxsplit=109)[-1] + '\t' + line_group3.split(maxsplit=109)[
        -1]

    fh_output.write(line_to_write + '\n')

    # Read in actual data line for the while loop to start
    line_group1 = fh_group1.readline().strip()
    line_group2 = fh_group2.readline().strip()
    line_group3 = fh_group3.readline().strip()

    count = 0   # For console output
    while line_group1 != '' and line_group2 != '' and line_group3 != '':
        # when end of file is not reached for all three files
        lst_g1 = line_group1.split()
        # Remove the first 109 columns by setting maxsplit in split()
        # Those columns are repeated (with 100 duplicated individuals)
        lst_g2 = line_group2.split(maxsplit=109)
        lst_g3 = line_group3.split(maxsplit=109)
        # Get snp ID of each group and compare
        snp_group1 = lst_g1[2]
        snp_group2 = lst_g2[2]
        snp_group3 = lst_g3[2]

        # Compare snp ids with the snp_keep list file (lst_snp_to_keep)
        # Pop snps from the list until it's empty
        if lst_snp_to_keep != []:
            snp_to_keep = lst_snp_to_keep.pop(0)
            while snp_group1 != snp_to_keep:
                # To see if current sno in the vcf.dose.gz file needs to be kept
                # Continue read the file until find the snp
                line_group1 = fh_group1.readline().strip()
                lst_g1 = line_group1.split()
                snp_group1 = lst_g1[2]
            while snp_group2 != snp_to_keep:
                line_group2 = fh_group2.readline().strip()
                lst_g2 = line_group2.split(maxsplit=109)
                snp_group2 = lst_g2[2]
            while snp_group3 != snp_to_keep:
                line_group3 = fh_group3.readline().strip()
                lst_g3 = line_group3.split(maxsplit=109)
                snp_group3 = lst_g3[2]

            # Write into output file
            # Piece together 3 files (first 100 individuals are replicated)
            line_to_write = line_group1 + '\t' + lst_g2[-1] + '\t' + lst_g3[-1]

            # print('group1:\t\tline len = ', len(lst_g1))
            # print('group2:\t\tline len = ', len(lst_g2[-1].split()))
            # print('group3:\t\tline len = ', len(lst_g3[-1].split()))
            # print(len(lst_g3))
            # quit()
            # print('line to write:\tlen =',len(line_to_write.split()), '\n-------------')

            fh_output.write(line_to_write + '\n')

        # Below is to keep console busy
        count = count + 1
        if count % 100 == 0:
            print('.', end='', flush=True)

        if count % 10000 == 0:
            print('\n' + chromosome + ':', count, 'SNPs have been merged', flush=True)
        elif count % 1000 == 0:
            print('', flush=True)

        line_group1 = fh_group1.readline().strip()
        line_group2 = fh_group2.readline().strip()
        line_group3 = fh_group3.readline().strip()

    print('\n' + chromosome + ':', count, 'SNPs have been merged')
    print('\n------------------------------------\nDone')

# ----------------------- End of helper functions -----------------------



