# Step 01
# This code is to check Rsq values across all input file .info.gz
# Create two files of variants to be excluded and kept
# Below information is saved in .txt files (variant_excluded.txt and variant_kept.txt):
#   - SNP, REF(0), ALT(1), and ALT_freq, MAF, Rsq (r2) per each group
#   - Each column is annotated with suffix "_group1", "_group2", etc.
#     (E.g. ALT_Freq_group1, MAF_group1, Rsq_group1)
# If other info is desired, change cols_to_keep parameter. Possible columns are:
#   - SNP, REF(0), ALT(1), ALT_Frq, MAF, AvgCall, Rsq, Genotyped
#   - LooRsq, EmpR, EmpRsq, Dose0, Dose1 (These seem to be NA for most samples)

import pandas as pd

# This function returns a list of DataFrames read from .info.gz files
# Parameter:
# - dict_flag: a dictionary containing flags and values. Such as '--missing':1, etc.
def __get_lst_info_df(dict_flags):
    # Get a list of names for .info.gz file based on names of input files (change suffix)
    # .info.gz and .dose.vcf.gz should be the same beside suffix
    print('\nBelow .info.gz files will be used:')
    lst_info_fn = []
    for fn in dict_flags['--input']:
        info_fn = fn.split('.')[0] + '.info.gz'
        lst_info_fn.append(info_fn)
        print('\t' + info_fn)

    lst_info_df = []  # A list to store .info.gz DataFrames
    for info_fn in lst_info_fn:
        try:
            df = pd.read_csv(info_fn, sep='\t', compression='gzip', dtype='str')
            lst_info_df.append(df)
        except:
            print('Error: File not found:', info_fn, '\n')
            raise IOError('File not found: ' + info_fn)
    return lst_info_df

# This function merge all .info.gz files (outer merge)
# Parameters:
# - lst_info_df: a list containing names of .info.gz files
# - cols_to_keep: a list of columns to be saved in the output file.
#                 In this version use ['SNP', 'REF(0)', 'ALT(1)', 'Genotyped', 'MAF', 'Rsq']
# Return: a Dataframe with all variants and columns REF(0), ALT(1), Genotyped, MAF and Rsq
# - For duplicated fields MAF and Rsq, column are renamed as MAF_group1, Rsq_group1, etc
def __merge_snps(lst_info_df, cols_to_keep=['SNP', 'REF(0)', 'ALT(1)', 'Genotyped', 'MAF', 'Rsq']):
    # REF(0), ALT(1), Genotyped should be the same for all files,
    # so only keep then in the first dataframe to merge
    cols_to_keep_v2 = cols_to_keep.copy()
    for val in ['REF(0)', 'ALT(1)', 'Genotyped']:
        cols_to_keep_v2.remove(val)

    i=0
    df_merged='' # A DataFrame of merged .info.gz df
    while i<len(lst_info_df):
        if i==0: # Merge the 1st and 2nd df
            df_merged = lst_info_df[i][cols_to_keep].merge(lst_info_df[i+1][cols_to_keep_v2],
                                                          how='outer', left_on='SNP', right_on='SNP',
                                                          suffixes=('_group'+str(i+1), '_group'+str(i+2)))
            i = i+2
        else:
            # Rename columns of df with correct suffix (ie. _groupX)
            df_merged = df_merged.merge(lst_info_df[i][cols_to_keep_v2].rename(columns={'MAF':'MAF_group'+str(i+1), 'Rsq':'Rsq_group'+str(i+1)}),
                                      how='outer', left_on='SNP',
                                      right_on='SNP')
            i+=1
    return df_merged
# ----------------- End of helper functions -----------------

# This function save excluded and kept variants into .txt files
# Parameter:
#  - df_merged: A DataFrame of variants from all input files with desired fields
#  - missing: A variant is allowed to be missing in this number of files. Otherwise it will be excluded
# Output:
#  - Two text files saved in current directory: variants_kept.txt and variants_excluded.txt
def __process_output(df_merged, dict_flags):
    # Set output file names
    to_keep_fn = 'variants_kept.txt'
    to_exclude_fn = 'variants_excluded.txt'
    missing = dict_flags['--missing']
    if missing==0: # Only save variants shared by all input files
        # Check columns of MAF (Rsq also works) such as: MAF_group1, MAF_group2,...
        # If MAF_groupX is NaN, then this variant is missing from the corresponding group.
        mask_to_keep = ''
        mask_to_exclude = ''
        for i in range(len(dict_flags['--input'])):
            col_name = 'MAF_group' + str(i + 1)
            if i == 0:  # The first iteration
                mask_to_keep = df_merged[col_name].notna()
                mask_to_exclude = df_merged[col_name].isna()
            else:
                mask_to_keep = mask_to_keep & df_merged[col_name].notna()
                mask_to_exclude = mask_to_exclude | df_merged[col_name].isna()
            df_merged[mask_to_keep].to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
            df_merged[mask_to_exclude].to_csv(to_exclude_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        print('\tNumber of saved variants:', df_merged[mask_to_keep].shape[0])
        print('\tNumber of excluded variants:', df_merged[mask_to_exclude].shape[0])
    else: # User specify allowed missingness, max value is number of input files
        lst_col_names = []
        for i in range(len(dict_flags['--input'])):
            lst_col_names.append('MAF_group' + str(i + 1))    # Save column names of MAF_group1, MAF_group2, etc.

        # Note thresh in dropna(): Require that many non-NA values in order to not drop
        inx_to_keep = df_merged[lst_col_names].dropna(thresh=len(dict_flags['--input'])-dict_flags['--missing']).index
        df_merged.iloc[inx_to_keep, :].to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        df_merged.drop(index=inx_to_keep).to_csv(to_exclude_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        print('\tNumber of saved variants:', len(inx_to_keep))
        print('\tNumber of excluded variants:', df_merged.drop(index=inx_to_keep).shape[0])

# A wrapper function to run this script
# Returns dict_flags for next step
def get_snp_list(dict_flags):
    lst_info_df = __get_lst_info_df(dict_flags)
    df_merged = __merge_snps(lst_info_df)
    print('\nNumber of variants:')
    print('\tTotal number from all input files:', df_merged.shape[0])
    __process_output(df_merged, dict_flags)
