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
import gzip

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

# This function merge all .info.gz files (outer merge) into a master Dataframe
# Parameters:
# - lst_info_df: a list containing names of .info.gz files
# - cols_to_keep: a list of columns to be saved in the output file.
#                 In this version use ['SNP', 'REF(0)', 'ALT(1)', 'Genotyped','ALT_Frq', 'MAF', 'Rsq']
# Return: a Dataframe with all variants and columns REF(0), ALT(1), Genotyped, ALT_Frq, MAF and Rsq
# - For duplicated fields MAF and Rsq, column are renamed as MAF_group1, Rsq_group1, etc
def __merge_snps(lst_info_df, cols_to_keep=['SNP', 'REF(0)', 'ALT(1)', 'Genotyped', 'ALT_Frq', 'MAF', 'Rsq']):
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
            df_merged = df_merged.merge(lst_info_df[i][cols_to_keep_v2].rename(columns={'MAF':'MAF_group'+str(i+1), 'Rsq':'Rsq_group'+str(i+1), 'ALT_Frq':'ALT_Frq_group'+str(i+1)}),
                                      how='outer', left_on='SNP',
                                      right_on='SNP')
            i+=1
    return df_merged

# This function calculated weighted r2 and MAF,
# assign values to a column (col_name_r2_combined, col_name_maf_combined) in df_merged
# Parameters:
# - df_merged: merged dataframe from all .info.gz files
# - col_name_r2_combined: column name of combined r2
# - col_name_maf_combined: column name of combined minor allele frequency
# - col_name_alt_frq: column name of combined alternative allele frequency (ALT_Frq)
# - lst_input_fn: a list of names of input files
def __calculate_weighted_r2_maf_altFrq(df_merged, col_name_r2_combined, col_name_maf_combined, col_name_alt_frq_combined, lst_input_fn):
    lst_number_of_individuals = []  # A list to store number of individuals of each input file
    # Count number of individuals in each input file
    for fn in lst_input_fn:
        with gzip.open(fn, 'rt') as fh:
            line = fh.readline().strip()
            while line[0:2]=='##':
                line = fh.readline().strip() # Read and skip header lines
            lst_tmp = line.split()
            start_pos = lst_tmp.index('FORMAT')
            total_len = len(lst_tmp)
            lst_number_of_individuals.append(total_len-(start_pos+1))

    # Calculate weighted r2 use this equation:
    # r2_combined = sum(r2_groupX * number_of_individuals_groupX) / sum(number_of_individuals_groupX)
    lst_col_names_r2_adj = []   # A list to store column names of r2 * weight
    lst_col_names_maf_adj = []  # For MAF, similar to lst_col_names_r2_adj
    lst_col_names_alt_frq_adj = [] # For ALT_Frq, same calculation as MAF

    # A list to store column names of r2 * weight/r2
    # In this way weight can be adjusted to reflect missing values in r2
    lst_col_names_weight_adj = []

    for i in range(len(lst_input_fn)):
        col_name_r2 = 'Rsq_group'+str(i+1)
        col_name_maf = 'MAF_group' + str(i + 1)
        col_name_alt_frq = 'ALT_Frq_group'+str(i+1)

        # !!! Must use pd.DataFrame.sum(), mean(), etc. functions to deal with missing values (NaN)
        # Otherwise will get NaN if adding NaN to another value directly.
        df_merged[col_name_r2+'_adj'] = df_merged[col_name_r2] * lst_number_of_individuals[i] # Rsq * # individuals
        df_merged[col_name_maf+'_adj'] = df_merged[col_name_maf] * lst_number_of_individuals[i]   # MAF * # individuals
        df_merged[col_name_alt_frq + '_adj'] = df_merged[col_name_alt_frq] * lst_number_of_individuals[i]  # MAF * # individuals
        # Weight is number of individuals
        df_merged[col_name_r2 + '_weight_adj'] = lst_number_of_individuals[i] * df_merged[col_name_r2]/ df_merged[col_name_r2]
        lst_col_names_r2_adj.append(col_name_r2+'_adj')
        lst_col_names_maf_adj.append(col_name_maf + '_adj')
        lst_col_names_alt_frq_adj.append(col_name_alt_frq + '_adj')
        lst_col_names_weight_adj.append(col_name_r2 + '_weight_adj')

    df_merged[col_name_alt_frq_combined] = df_merged[lst_col_names_alt_frq_adj].sum(axis=1) / sum(lst_number_of_individuals)
    df_merged[col_name_maf_combined] = df_merged[lst_col_names_maf_adj].sum(axis=1) / sum(lst_number_of_individuals)
    df_merged[col_name_r2_combined] = df_merged[lst_col_names_r2_adj].sum(axis=1) / df_merged[lst_col_names_weight_adj].sum(axis=1)

    # Drop columns that are not needed in output
    df_merged.drop(labels=lst_col_names_r2_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_maf_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_alt_frq_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_weight_adj, axis=1, inplace=True)

    return lst_number_of_individuals    # Return number of individuals in each file
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
    mask_to_keep = ''   # Use this when --missing is 0
    mask_to_exclude = ''    # Use this when --missing is 0

    lst_col_names = []  # A list to store column names of Rsq (such as Rsq_group1, Rsq_group2, etc.)

    for i in range(len(dict_flags['--input'])):
        col_name_r2 = 'Rsq_group' + str(i + 1)
        col_name_maf = 'MAF_group' + str(i + 1)
        col_name_alt_frq = 'ALT_Frq_group' + str(i + 1) # Column name of alternative allele
        # Convert values of r2, MAF and ALT_Frq columns from strings to numbers for later calculations
        df_merged[col_name_r2] = pd.to_numeric(df_merged[col_name_r2], errors='coerce')
        df_merged[col_name_maf] = pd.to_numeric(df_merged[col_name_maf], errors='coerce')
        df_merged[col_name_alt_frq] = pd.to_numeric(df_merged[col_name_alt_frq], errors='coerce')

        if missing == 0:  # Only save variants shared by all input files
            # Check columns of Rsq (MAF also works) such as: Rsq_group1, Rsq_group2,...
            # If MAF_groupX is NaN, then this variant is missing from the corresponding group.
            if i == 0:  # The first iteration
                mask_to_keep = df_merged[col_name_r2].notna()
                mask_to_exclude = df_merged[col_name_r2].isna()
            else:
                mask_to_keep = mask_to_keep & df_merged[col_name_r2].notna()
                mask_to_exclude = mask_to_exclude | df_merged[col_name_r2].isna()

            # Check if need to filter for SNPs using given r2 threshold
            if dict_flags['--r2_threshold'] != 0:
                mask_to_keep = mask_to_keep & (df_merged[col_name_r2] >= dict_flags['--r2_threshold'])
                mask_to_exclude = mask_to_exclude | (df_merged[col_name_r2] < dict_flags['--r2_threshold'])

        lst_col_names.append(col_name_r2)  # Save column names of MAF_group1, MAF_group2, etc.

    # Calculate combined r2 (first, weighted_average or mean)
    lst_number_of_individuals = []  # Print this info if '--r2_output' is 'weighted_average':
    col_name_r2_combined = 'Rsq_combined'
    col_name_maf_combined = 'MAF_combined'
    col_name_alt_frq_combined = 'ALT_Frq_combined'
    if dict_flags['--r2_output'] == 'first':  # use r2 of the first input file
        df_merged[col_name_r2_combined] = df_merged['Rsq_group1']
    elif dict_flags['--r2_output'] == 'weighted_average':
        lst_number_of_individuals = __calculate_weighted_r2_maf_altFrq(df_merged, col_name_r2_combined,
                                                                       col_name_maf_combined, col_name_alt_frq_combined,
                                                                       dict_flags['--input'])
    else:  # Mean
        df_merged[col_name_r2_combined] = df_merged[lst_col_names].mean(axis=1)

    # Save result to output files
    if missing == 0:
        df_merged[mask_to_keep].to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        df_merged[mask_to_exclude].to_csv(to_exclude_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        print('\tNumber of saved variants:', df_merged[mask_to_keep].shape[0])
        print('\tNumber of excluded variants:', df_merged[mask_to_exclude].shape[0])
    else:   # Use user specify allowed missingness, max value is number of input files
        # Note thresh in dropna(): Require that many non-NA values in order to not drop
        inx_to_keep = df_merged[lst_col_names].dropna(
            thresh=len(dict_flags['--input']) - dict_flags['--missing']).index
        df_merged.iloc[inx_to_keep, :].to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'])
        df_merged.drop(index=inx_to_keep).to_csv(to_exclude_fn, index=False, sep='\t',
                                                 na_rep=dict_flags['--na_rep'])
        print('\tNumber of saved variants:', len(inx_to_keep))
        print('\tNumber of excluded variants:', df_merged.drop(index=inx_to_keep).shape[0])

    print('\nNumbers of individuals in each input file:', lst_number_of_individuals)

# A wrapper function to run this script
# Returns dict_flags for next step
def get_snp_list(dict_flags):
    lst_info_df = __get_lst_info_df(dict_flags)
    df_merged = __merge_snps(lst_info_df)
    print('\nNumber of variants:')
    print('\tTotal number from all input files:', df_merged.shape[0])
    __process_output(df_merged, dict_flags)

