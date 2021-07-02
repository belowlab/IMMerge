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



# --------------------------------------------
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












# Close all file handles when done
for fh in lst_input_fh: fh.close()
fh_output.close()