# This code is used to rename vairant ID of merged files
# TopMED imputation files name SNPs as chr9:12065:C:A
# Need to rename to format as 'rs00000' for HLA analysis


import os
os.chdir('/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/')

def get_snp_and_info_lst(
        variants_kept_fn='/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/variants_kept.txt'):
    df_snps = pd.read_csv(variants_kept_fn, sep='\t')
    lst_snp = df_snps['SNP']
    lst_alt_frq = df_snps['ALT_Frq_combined']
    lst_maf = df_snps['MAF_combined']
    lst_r2 = df_snps['Rsq_combined']
    return (lst_snp, lst_alt_frq, lst_maf, lst_r2)


import gzip
import pandas as pd

print('\nStart merging files:')
number_of_dup_id = 10  # User provided number of duplicated IDs in each sub input file
lst_input_fn = [
    '/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group1.dose.vcf.gz',
    '/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group2.dose.vcf.gz',
    '/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/sample_data/subset_short_chr21_group3.dose.vcf.gz']
output_fn = 'chr21_test.dose.vcf.gz'
fh_output = gzip.open(output_fn, 'wt')  # Open for writing (appending)

lst_input_fh = []  # A list to store file handles of input files
for fn in lst_input_fn: lst_input_fh.append(gzip.open(fn, 'rt'))

# Read in SNPs (with corresponding combined MAF and r2 )need to be kept from file
lst_snp, lst_alt_frq, lst_maf, lst_r2 = get_snp_and_info_lst()

# ---------------------DONE!!!!  Write header lines-------------------------
# Read trough header lines (line start with ##) and write into output file
# Process all files at the same time. All input files should have the same number of header lines
line = ''
for fh in lst_input_fh:  # Read the 1st line of all input files
    line = fh.readline().strip()
fh_output.write(line)
line = lst_input_fh[0].readline().strip()  # Read the 2nd line of the first input files

inx_indiv_id_starts = 0  # From which column genotype data starts
inx_snp_ID = 0  # Column index of SNP id (in current TOPMed version is 2(2021/0708))
while line[0] == '#':
    if line[0:2] == '##':  # Info lines
        fh_output.write('\n' + line)
        for fh in lst_input_fh[1:]:  # Read the rest fh
            line = fh.readline().strip()
    else:  # Merge and write header line (line starts with single #)
        # In current TOPMed version these are columns of shared information, then followed by individual IDs:
        #   - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
        inx_indiv_id_starts = line.split().index('FORMAT') + 1  # Find from which column genotype data starts
        inx_snp_ID = line.split().index('ID')   # Find column index of SNP ID (such as ID chr21:10000931:T:G)
        print('\tIndividual genotype data starts from column #:', inx_indiv_id_starts + 1)

        # Create merged header line, initialize with header from the first file
        line_to_write = '\n' + line
        for fh in lst_input_fh[1:]:
            line = fh.readline().strip()
            # Set maxsplit in slit() function to remove duplicated IDs
            line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
        fh_output.write(line_to_write + '\n')
        break  # Break at the column header line
    line = lst_input_fh[0].readline().strip()  # Read line of the first file to decide what to do
# -------------------------- DONE!!!! write header lines -----------------------------


# -------------------- Remove below in real code ------------------------
inx_indiv_id_starts = 9
inx_snp_ID = 2
inx_alt_frq_col = 4
inx_info_column = 7
inx_snp_id_column = 2
lst_number_of_individuals = [21, 21, 21]
dict_flags = dict()
dict_flags['--input'] = ['sample_data/subset_short_chr21_group1.dose.vcf.gz',
                         'sample_data/subset_short_chr21_group2.dose.vcf.gz',
                         'sample_data/subset_short_chr21_group3.dose.vcf.gz']
dict_flags['--missing'] = 1
dict_flags['--r2_threshold'] = 0.2
dict_flags['--r2_output'] = 'weighted_average'
dict_flags['--duplicate_id'] = 10
dict_flags['--output'] = 'merged'
dict_flags['--verbose'] = 'true'
dict_flags['--na_rep'] = 'NA'


# -------------------- Remove above in real code ------------------------


# -------------------------------- Merge rest lines-----------------------------------

# This function changes INFO column using combined MAF, alt allele freq and r2
# INFO value format is: AF=xxx;MAF=xxx;R2=xxx;IMPUTED/TYPED/TYPED_ONLY
# Return updated value of INFO of a single SNP
def change_INFO_val(org_INFO, alt_frq_combined, MAF_combined, r2_combined):
    genotype_type = org_INFO.split(';')[-1]
    updated_INFO = 'AF=' + alt_frq_combined + ';MAF=' + MAF_combined + ';R2=' + r2_combined + ';' + genotype_type
    return updated_INFO

# This function read line of a given file handle, until find line of given SNP
# Return string of that line
def search_SNP_and_read_lines(snp, fh):
    line = fh.readline().strip()  # Read in each line to search for the given variant
    input_snp = line.split(maxsplit=inx_snp_ID + 1)[-2]
    while input_snp != snp and line != '': # Keep read line until find given snp
        input_snp = line.split(maxsplit=inx_snp_ID + 1)[-2]
        line = fh.readline().strip()  # Read in each line to search for the given variant
    return line

# This function takes in a variant ID,
# check every input file to find corresponding line and return a merged line.
# If a snp is not found in certain input file, fill with symbol of missing values for that file
# Missing values symbol can be any user supplied string, default is 'NA',
def merge_individual_variant(snp, lst_alt_frq_val, lst_maf_combined, lst_r2_combined, lst_input_fh):
    merged_line = ''
    # flag_first_na: An indicator to show whether a variant is missing from the first input file
    # If missing, need to fill info columns from other input file
    flag_first_na = False
    for i in range(len(lst_input_fh)):
        if i==0:    # If the first input file
            if lst_alt_frq_val[i] != dict_flags['--na_rep']:
                # If this is the first input file and variant is not missing,
                # then keep all columns including info columns
                merged_line = search_SNP_and_read_lines(snp, lst_input_fh[i])
            else:
                # If variant is missing from the file, then fill with 'NA' and exclude the first few info columns
                # Need to get values of info columns from other input files
                merged_line = '\t'.join([dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i])])
                flag_first_na = True # Indicate flag: Fill info columns later with other input files
        else: # If this is not the first file
            if lst_alt_frq_val[i] == dict_flags['--na_rep']:
                # If variant is missing from this input file, then fill 'NA' or user provided missing value symbol
                merged_line = merged_line + '\t' + '\t'.join([dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i])])
            else: # If not missing, split and remove info columns
                line = search_SNP_and_read_lines(snp, lst_input_fh[i])
                merged_line = merged_line + '\t' + line.split(maxsplit=inx_indiv_id_starts)[-1]
                if flag_first_na:   # If info columns is still missing, add these columns back
                    info_cols = '\t'.join(line.split(maxsplit=inx_indiv_id_starts)[:-1])
                    merged_line = info_cols + '\t' + merged_line
                    flag_first_na = False
    return merged_line


count = 0  # For console output
# Load variants_kept.txt as a reference of merged file
# Compare snp ids in each input file with the SNPs in variants_kept.txt
fh_snp_kept = open('variants_kept.txt')
fh_snp_kept.readline()  # Skip the first line (header)

# Remove the first 9 columns and duplicate_id columns by setting maxsplit in split()
while True:
    line_snp_kept = fh_snp_kept.readline().strip()
    if line_snp_kept != '':
        lst_snp_kept = line_snp_kept.split()
        snp = lst_snp_kept[0]  # Get ID of SNP to be kept
        # Get ALT_frq, MAF and r2 values
        # Use ALT_frq_groupX values to decide whether SNP is missing in the Xth file
        lst_alt_frq_val = []
        lst_maf_combined = []
        lst_r2_combined = []

        inx_alt_frq_col = 4    # Collect index of ALT_frq column
        inx_maf_col = 5  # Collect index of ALT_frq column
        inx_r2_col = 6  # Collect index of ALT_frq column
        for i in range(len(lst_input_fh)):
            lst_alt_frq_val.append(lst_snp_kept[inx_alt_frq_col])
            lst_maf_combined.append(lst_snp_kept[inx_maf_col])
            lst_r2_combined.append(lst_snp_kept[inx_r2_col])
            inx_alt_frq_col += 3
            inx_maf_col += 3
            inx_r2_col += 3

        # Search line of this snp in each input file, fill with missing values if not exist
        # Use ALT_Frq_groupX column in variant_kept.txt file to indicate if this snp exist in input files
        # Eg. if ALT_Frq_group1 is missing, then this snp does not exist in the first input file
        print('SNP to be kept:', snp)
        merged_line = merge_individual_variant(snp, lst_alt_frq_val, lst_maf_combined, lst_r2_combined, lst_input_fh)
        count += 1
        # print('.', end='', flush=True)
        # print('{} variants merged'.format(count))
    else:
        print('End of snp to be kept file')
        break  # If reach end of SNP to be kept, then stop

    print(lst_alt_frq_val)
    print(merged_line)
    print('====================')


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
