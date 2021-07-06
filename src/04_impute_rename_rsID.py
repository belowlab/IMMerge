# This code is used to rename vairant ID of merged files
# TopMED imputation files name SNPs as chr9:12065:C:A
# Need to rename to format as 'rs00000' for HLA analysis

# 2020/09/18 For now only need to rename SNPs on chr6 (HLA type) and chr9 (ABO blood type)
def get_snp_and_info_lst(variants_kept_fn='/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/variants_kept.txt'):
    df_snps = pd.read_csv(variants_kept_fn, sep='\t')
    lst_snp = df_snps['SNP']
    lst_alt_frq = df_snps['ALT_Frq_combined']
    lst_maf = df_snps['MAF_combined']
    lst_r2 = df_snps['Rsq_combined']
    return (lst_snp, lst_alt_frq, lst_maf, lst_r2)


import gzip
import pandas as pd
print('\nStart merging files:')
number_of_dup_id = 10 # User provided number of duplicated IDs in each sub input file
lst_input_fn = ['/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group1.dose.vcf.gz',
                '/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group2.dose.vcf.gz',
                '/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group3.dose.vcf.gz']
output_fn = 'chr21_test.dose.vcf.gz'
fh_output = gzip.open(output_fn, 'wt')  # Open for writing (appending)

lst_input_fh = []  # A list to store file handles of input files
for fn in lst_input_fn: lst_input_fh.append(gzip.open(fn, 'rt'))

# Read in SNPs (with corresponding combined MAF and r2 )need to be kept from file
lst_snp, lst_alt_frq, lst_maf, lst_r2 = get_snp_and_info_lst()




# ---------------------DONE!!!!-------------------------
# Read trough header lines (line start with ##) and write into output file
# Process all files at the same time. All input files should have the same number of header lines
line = ''
for fh in lst_input_fh: # Read the 1st line of all input files
    line = fh.readline().strip()
fh_output.write(line)
line = lst_input_fh[0].readline().strip()   # Read the 2nd line of the first input files

inx_indiv_id_starts = 0 # From which column genotype data starts
while line[0] == '#':
    if line[0:2] == '##': # Info lines
        fh_output.write('\n'+line)
        # line = fh.readline().strip()
        # if line[0:2] == '##':  # Info lines
        for fh in lst_input_fh[1:]: # Read the rest fh
            line = fh.readline().strip()
    else: # Merge and write header line (line starts with single #)
        # In current TOPMed version these are columns of shared information, then followed by individual IDs:
        #   - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
        inx_indiv_id_starts = line.split().index('FORMAT')+1 # Find from which column genotype data starts
        print('\tIndividual genotype data starts from column #:', inx_indiv_id_starts + 1)
        # Create merged header line, initialize with header from the first file
        line_to_write = '\n' + line
        for fh in lst_input_fh[1:]:
            line = fh.readline().strip()
            # Set maxsplit in slit() function to remove duplicated IDs
            line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
        fh_output.write(line_to_write + '\n')
        break

    line = lst_input_fh[0].readline().strip()   # Read line of the first file to decide what to do





# -------------------- Remove these in real code ------------------------
inx_indiv_id_starts = 9
# inx_flag_col = 4
inx_info_column = 7
inx_snp_id_column = 2

dict_flags = dict()
dict_flags['--input'] = ['sample_data/subset_short_chr21_group1.dose.vcf.gz',
                         'sample_data/subset_short_chr21_group2.dose.vcf.gz',
                         'sample_data/subset_short_chr21_group3.dose.vcf.gz']
dict_flags['--missing'] = 1
dict_flags['--r2_threshold'] = 0.2
dict_flags['--r2_output'] = 'weighted_average'
dict_flags['--duplicate_id'] = 10
dict_flags['--output']= 'merged'
dict_flags['--verbose'] = 'true'
dict_flags['--na_rep'] = 'NA'

# -------------------- Remove this in real file ------------------------


# Merge rest lines------------------------------------------
count = 0   # For console output
# Load variants_kept.txt as a reference of merged file
# Compare snp ids in each input file with the SNPs in variants_kept.txt
fh_snp_kept = open('variants_kept.txt')
fh_snp_kept.readline()  # Skip the first line (header)

# Remove the first 9+duplicate_id columns by setting maxsplit in split()
while True:
    line_snp_kept = fh_snp_kept.readline().strip()
    if line_snp_kept != '':
        snp = line_snp_kept.split()[0]  # Get SNP ID for comparison
        print('SNP to be kept:', snp)
    else:
        print('End of snp to be kept file')
        break # If reach end of SNP to be kept, then stop

    # Replace alt_frq, MAF and r2 in INFO column with averaged values

    # # Use value of AlT_frq_groupX to decide how to read each input file
    # inx_flag_col = 4 # Column index of AlT_frq_group1
    # for i in range(len(lst_input_fh)): # check SNP in each input file
    #     line = lst_input_fh[i].readline().strip()
    #     if line != '': # has not reach end of current input file
    #         current_snp = line.split()[inx_snp_id_column]
    #         line_to_write = ''
    #         missing_info_cols = False # In case the first input file misses this SNP, need to get info cols from other input files
    #
    #         if line_snp_kept.split()[inx_flag_col] != dict_flags['--na_rep']:
    #             # If value at AlT_frq_groupX is not NA (missing),
    #             # then keep reading current input file until the SNP is found
    #             if current_snp == snp:
    #                 if i==0: # If the first input file, need to retain info columns
    #                     line_to_write = '\n' + line # Output merged line
    #                 else:
    #                     line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
    #             else:   # Keep reading until find the SNP
    #                 while current_snp != snp:
    #                     line = lst_input_fh[i].readline().strip()
    #                     current_snp = line.split()[inx_snp_id_column]
    #                 if missing_info_cols: # Add missing info columns if missing
    #                     info_cols = '\t'.join(line.split(maxsplit=inx_indiv_id_starts)[:inx_indiv_id_starts])
    #                     line_to_write = '\n' + info_cols + '\t' + line_to_write
    #                     missing_info_cols = False
    #         else:
    #         # If value at AlT_frq_groupX is NA, then skip reading current input file,
    #         # since the SNP is missing from this file
    #         # Append NA (missing) to the line
    #         if i == 0:  # If the SNP is missing from the first file, then need to use info column from other files
    #             # Save spot for missing info
    #             line_to_write = '\t'.join([dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i])])
    #             missing_info_cols = True
    #         else:
    #             if missing_info_cols: # If info columns are still missing, get them from current file
    #                 info_cols = '\t'.join(line.split(maxsplit=inx_indiv_id_starts)[:inx_indiv_id_starts])
    #                 line_to_write = '\n'+info_cols + '\t' + line_to_write
    #                 missing_info_cols = False
    #             line_to_write = line_to_write + '\t' + '\t'.join([dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i])])
    #
    #     fh_output.write(line_to_write)






# Close all file handles when done
for fh in lst_input_fh: fh.close()
fh_output.close()