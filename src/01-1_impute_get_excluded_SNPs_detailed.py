# This code is to check Rsq values accross all three groups of imputation file .info.gz
# Want to do this in SNP excluded file, output below info:
#   - SNP, REF(0), ALT(1), and ALT_freq, MAF, Rsq per each group
#       (E.g. ALT_Freq_group1, MAF_group1, Rsq_group1)
# This is a similar code to 01_impute_get_keep_or_exclude_snp_list.py.
# But only look for exclude SNPs and output a detailed file

import pandas as pd

# ----------------------- Helper functions -----------------------
# This function reads in a file with one snp in each line, (duplication of the same function in 02_impute_merge_files.py)
# then return a list of all SNPs.
# For this version, the first header line 'SNP' need to be skipped
def get_snp_lst(snp_fn):
    lst_snp = []
    with open(snp_fn, 'r') as fh:
        line = fh.readline()  # Skip the header line
        line = fh.readline().strip()
        while line != '':
            lst_snp.append(line)
            line = fh.readline().strip()
    return lst_snp

def get_excluded_snps(chromosome = 'chr1', info_gz_dir = '/data100t1/share/BioVU/TOPMed_imputation_072020/',
                      output_dir='/data100t1/home/wanying/imputation/snp_keep_or_exclude_list/'):
    # Create a list of .info.gz file names
    chr_info_fn_lst = [chromosome + '_group1.info.gz', chromosome + '_group2.info.gz', chromosome + '_group3.info.gz']

    df_group1 = pd.read_csv(info_gz_dir + chr_info_fn_lst[0], sep='\t', compression='gzip', dtype='str')
    df_group2 = pd.read_csv(info_gz_dir + chr_info_fn_lst[1], sep='\t', compression='gzip', dtype='str')
    df_group3 = pd.read_csv(info_gz_dir + chr_info_fn_lst[2], sep='\t', compression='gzip', dtype='str')

    df_merge = df_group1.iloc[:, [0, 1, 2, 3, 4, 6]].merge(df_group2.iloc[:, [0, 3, 4, 6]], how='outer',
                                                           left_on='SNP', right_on='SNP', indicator=True,
                                                           suffixes=('_group1', '_group2'))
    df_merge.rename(columns={'_merge': 'merge_group1_group2'}, inplace=True)
    df_merge = df_merge.merge(df_group3.iloc[:, [0, 3, 4, 6]], left_on='SNP', right_on='SNP', how='outer',
                              indicator=True)
    df_merge.rename(columns={'ALT_Frq': 'ALT_Frq_group3',
                             'MAF': 'MAF_group3',
                             'Rsq': 'Rsq_group3',
                             '_merge': 'merge_g1g2_group3'},
                    inplace=True)

    mask_exclude = (df_merge['merge_group1_group2'] != 'both') | (df_merge['merge_g1g2_group3'] != 'both')
    df_merge = df_merge.loc[mask_exclude].drop(['merge_group1_group2', 'merge_g1g2_group3'], axis=1)

    output_exclude_fn = chromosome + '_merge_excluded_snps_detailed.txt'
    # output_keep_fn = chromosome + '_merge_keep_snps.txt'
    df_merge.to_csv(output_dir+output_exclude_fn, sep='\t', index=False, na_rep='NA')

# ----------------------- End of helper functions -----------------------

# get_excluded_snps(chromosome = 'chr21')

# chr_number_list = ['1','2','3','4','5','6','7','13','14','15','16','17','18','19','20', '21','22','X']
# chr_number_list = ['8', '9']
#
# for i in chr_number_list:
#     get_excluded_snps(chromosome='chr'+i)
get_excluded_snps(chromosome='chr11')