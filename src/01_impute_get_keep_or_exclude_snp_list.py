import pandas as pd
# ----------------------- Helper functions -----------------------

# This function scan trough .info.gz files (Does not need to unzip .info.gz files manually)
# Then save results in chr*_merge_excluded_snps.txt and chr*_merge_keep_snps.txt
# Parameters
#   - chr: chromosome to be processed, such as "chr1"
#   - input_dir, output_dir: directory of input and output files
# Output: (two output files)
#   - Store SNPs shared in all three groups in chr*_merge_keep_snps.txt
#   - Store other SNPs in chr*_merge_excluded_snps.txt. They will be removed in the final merged .dose.vcf.gz file
def get_excluded_snps(chromosome='chr1',
                      input_dir='/data100t1/share/BioVU/TOPMed_imputation_072020/',
                      output_dir='/data100t1/home/wanying/imputation/snp_keep_or_exclude_list/'):
    # Create a list of .info.gz file names
    chr_info_fn = [chromosome + '_group1.info.gz', chromosome + '_group2.info.gz', chromosome + '_group3.info.gz']

    # df = pd.read_csv(dir+chr6_info_fn[0], sep='\t', skiprows=18, compression='gzip')

    # Read in SNP info, these files are much smaller than .dose.vcf.gz files
    df_group1 = pd.read_csv(input_dir + chr_info_fn[0], sep='\t', compression='gzip')
    df_group2 = pd.read_csv(input_dir + chr_info_fn[1], sep='\t', compression='gzip')
    df_group3 = pd.read_csv(input_dir + chr_info_fn[2], sep='\t', compression='gzip')

    # Merge group1 and group2
    df_merge = df_group1[['SNP', 'MAF']].merge(df_group2[['SNP', 'MAF']], how='outer',
                                               left_on='SNP', right_on='SNP', indicator=True)
    # Drop unnecessary columns
    df_merge.drop(['MAF_x', 'MAF_y'], axis=1, inplace=True)
    df_merge.columns = ['SNP', 'merge_g1_g2']
    # Change merged indicator values, since there are 3 groups. Need another merge
    df_merge.replace('left_only', 'group1', inplace=True)
    df_merge.replace('right_only', 'group2', inplace=True)
    df_merge.replace('both', 'group1_group2', inplace=True)

    # Merge and drop unnecessary columns again
    df_merge = df_merge.merge(df_group3[['SNP', 'MAF']], how='outer', left_on='SNP',
                              right_on='SNP', indicator=True)
    df_merge.drop(['MAF'], axis=1, inplace=True)
    df_merge.columns = ['SNP', 'merge_g1_g2', 'merge_g1g2_g3']
    df_merge.replace('left_only', 'group1_group2', inplace=True)
    df_merge.replace('right_only', 'group3', inplace=True)
    df_merge.replace('both', 'all', inplace=True)

    # Prepare for output
    mask_exclude = (df_merge['merge_g1_g2'] != 'group1_group2') | (df_merge['merge_g1g2_g3'] != 'all')
    mask_keep = (df_merge['merge_g1_g2'] == 'group1_group2') & (df_merge['merge_g1g2_g3'] == 'all')

    output_exclude_fn = chromosome + '_merge_excluded_snps.txt'
    output_keep_fn = chromosome + '_merge_keep_snps.txt'
    # df_merge.loc[mask_exclude].drop(['merge_g1_g2', 'merge_g1g2_g3'], axis=1).to_csv(output_dir+'archive_excluded_SNPs_no_detail/' + output_exclude_fn,
    #                                                                                  sep='\t', index=False)
    df_merge.loc[mask_keep].drop(['merge_g1_g2', 'merge_g1g2_g3'], axis=1).to_csv(output_dir + output_keep_fn, sep='\t',
                                                                                  index=False)
# ----------------------- End of helper functions -----------------------


# chr_number_list = ['1','2','3','4','5','6','7','13','14','15','16','17','18','19','20','22','X']
#
# for i in chr_number_list:
#     get_excluded_snps(chromosome='chr'+i)

get_excluded_snps(chromosome='chr11')