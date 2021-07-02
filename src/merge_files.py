# Step 02
# This code is used to merge imputed files from TopMed based on list of saved variants in prior steps
# The goal is to combine post-imputation vcf.gz files into one vcf.gz file

print('\nVersion 1.0')

import pandas as pd
import gzip
import sys
import process_args
import get_SNP_list
import check_r2_setting_for_imputation
import time

start_time = time.time()    # Track execution time

# File name can be passed to this code in terminal, or use import this code as in a script (need to modify a little)
args = sys.argv
verbose=True
# Process terminal input
# dict_flags contains values for below flags:
#   --input, --output, --verbose, --missing, --r2_threshold, --r2_output
dict_flags = process_args.process_args(args) # Process terminal input
if dict_flags['--verbose']!='true': verbose=False

check_r2_setting_for_imputation.check_imputatilson_parameters(lst_fn=dict_flags['--input'])
lst_number_of_individuals = get_SNP_list.get_snp_list(dict_flags)


# ----------------------- Helper functions -----------------------
'''
# This function reads in a file with one snp in each line (ie. variants_kept.txt),
# then return 3 lists of all SNPs, ALT_Frq_combined, MAF_combined and Rsq_combined
# The first line needs to be skipped (header line)
def get_snp_and_info_lst(variants_kept_fn='variants_kept.txt'):
    df_snps = pd.read_csv(variants_kept_fn, sep='\t')
    lst_snp = df_snps['SNP']
    lst_alt_frq = df_snps['ALT_Frq_combined']
    lst_maf = df_snps['MAF_combined']
    lst_r2 = df_snps['Rsq_combined']
    return (lst_snp, lst_alt_frq, lst_maf, lst_r2)
'''


# Merge input files (.dose.vcf.gz) together
# Parameters:
# - lst_number_of_individuals: number of individuals in each input file
def merge_files(dict_flags=dict_flags, lst_number_of_individuals=lst_number_of_individuals):
    print('\nStart merging files:')
    number_of_dup_id = dict_flags['--duplicate_id'] # User provided number of duplicated IDs in each sub input file
    lst_input_fn = dict_flags['--input']
    output_fn = dict_flags['--output'] + '.dose.vcf.gz'
    fh_output = gzip.open(output_fn, 'wt')  # Open for writing (appending)

    lst_input_fh = []  # A list to store file handles of input files
    for fn in lst_input_fn: lst_input_fh.append(gzip.open(fn, 'rt'))

    # # Read in SNPs (with corresponding combined MAF and r2 )need to be kept from file
    # lst_snp, lst_alt_frq, lst_maf, lst_r2 = get_snp_and_info_lst()

    # Read trough header lines (line start with ##) and write into output file
    # Process all files at the same time. All input files should have the same number of header lines
    line = ''
    for fh in lst_input_fh:  # Read the 1st line of all input files
        line = fh.readline().strip()
    fh_output.write(line)
    line = lst_input_fh[0].readline().strip()  # Read the 2nd line of the first input files

    inx_indiv_id_starts = 0  # From which column genotype data starts
    inx_info_column = 0 # Find index of column "INFO", inorder to replace values later
    inx_snp_id_column = 0   # Find index of column "ID" (ie. column of SNP IDs)
    while line[0] == '#':
        if line[0:2] == '##':  # Info lines
            fh_output.write('\n' + line)
            # line = fh.readline().strip()
            # if line[0:2] == '##':  # Info lines
            for fh in lst_input_fh[1:]:  # Read the rest fh
                line = fh.readline().strip()
        else:  # Merge and write header line (line starts with single #)
            # In current TOPMed version these are columns of shared information, then followed by individual IDs:
            #   - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
            inx_indiv_id_starts = line.split().index('FORMAT') + 1  # Find from which column genotype data starts
            inx_info_column = line.split().index('INFO')
            inx_snp_id_column = line.split().index('ID')
            print('\tIndividual genotype data starts from column #:', inx_indiv_id_starts + 1)
            # Create merged header line, initialize with header from the first file
            line_to_write = '\n' + line
            for fh in lst_input_fh[1:]:
                line = fh.readline().strip()
                # Set maxsplit in slit() function to remove duplicated IDs
                line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
            fh_output.write(line_to_write + '\n')
            break

        line = lst_input_fh[0].readline().strip()  # Read line of the first file to decide what to do


    # Merge rest lines------------------------------------------
    count = 0   # For console output

    # Load variants_kept.txt as a reference of merged file
    # Compare snp ids in each input file with the SNPs in variants_kept.txt
    fh_snp_kept = open('variants_kept.txt')
    fh_snp_kept.readline()  # Skip the first line (header)

    # Remove the first 9+duplicate_id columns by setting maxsplit in split()

    while True:
        line_snp_kept = fh_snp_kept.readline().strip()
        snp = line_snp_kept.split()[0]  # Get SNP ID for comparison

        # line_to_write = '\n' + line # Output merged line
        lst_snp_id = line.split(maxsplit=inx_snp_id_column+1)[inx_snp_id_column]
        # Replace alt_frq, MAF and r2 in INFO column with averaged values

        inx_flag_col = 4 # Column index of AlT_frq_group1
        for i in range(len(lst_input_fh)): # check SNP in each input file
            line = lst_input_fh[i].readline().strip()
            current_snp = line.split()[inx_snp_id_column]
            if line_snp_kept.split()[inx_flag_col] != dict_flags['--na_rep']:
                # If value at AlT_frq_groupX is not NA (missing),
                # then keep reading current input file until the SNP is found
                if current_snp == snp:
                    if i==0: # If the first input file, need to retain info columns
                        line_to_write = '\n' + line # Output merged line
                    else:
                        line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
                else:   # Keep reading until find the SNP
                    while current_snp != snp:
                        line = lst_input_fh[i].readline().strip()
                        current_snp = line.split()[inx_snp_id_column]
            else:
                # If value at AlT_frq_groupX is NA, then skip reading current input file,
                # since the SNP is missing from this file
                # Append NA (missing) to the line
                if i == 0:  # If the SNP is missing from the first file, then need to use info column from other files
                    pass #!!!!!!!
                else:
                    number_of_individuals = 0
                    line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]

            print('\n',line.split()[inx_info_column].split(';'))
            print('\n2.',line_snp_kept)
            print('\n3.', line[:50])
            print('\n4.', current_snp, lst_input_fh[i])
            print('\n5.',line_snp_kept.split()[inx_flag_col], line_snp_kept.split()[inx_flag_col]=='2e-05')
            # line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]


            break
        # print(lst_snp_id)
        break
    # ------------------------------------------

    # Close all file handles when done
    for fh in lst_input_fh: fh.close()
    fh_output.close()
    fh_snp_kept.close()


merge_files()
duration = time.time() - start_time
hours = duration // 3600
duration = duration % 3600
minutes = duration // 60
seconds = duration % 60
print('\nDuration: {:.2f} hours, {:.2f} minutes, {:.2f} seconds'.format(hours, minutes, seconds))
print('End of run')


'''
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
'''
# ----------------------- End of helper functions -----------------------



