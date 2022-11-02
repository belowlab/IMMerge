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
import numpy as np
import os

# This function returns a list of DataFrames read from .info.gz files
# Parameter:
# - dict_flag: a dictionary containing flags and values. Such as '--missing':1, etc.
def __get_lst_info_df(dict_flags):
    # Get a list of names for .info.gz file based on names of input files (change suffix)
    # .info.gz and .dose.vcf.gz should be the same beside suffix
    print('\nBelow .info.gz files will be used:')
    for fn in dict_flags['--info']:
        print('\t' + fn)

    lst_info_df = []  # A list to store .info.gz DataFrames
    for info_fn in dict_flags['--info']:
        try:
            # Read in as string data type to avoid low memory warning (convert to numeric later)
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
#                 In this version use ['SNP', 'REF(0)', 'ALT(1)', 'Genotyped','ALT_Frq', 'MAF', 'Rsq'], last 3 items must in order AF, MAF, Rsq
#                 Info files generated from make_info.py follow this naming style
# Return: a Dataframe with all variants and columns REF(0), ALT(1), Genotyped, ALT_Frq, MAF and Rsq
# - For duplicated fields MAF and Rsq, column are renamed as MAF_group1, Rsq_group1, etc
# - lst_index_col_names: a list to store index column names of each dataframe, such as [index_group1, index_group2, ...]
def __merge_snps(lst_info_df, cols_to_keep=['SNP', 'REF(0)', 'ALT(1)', 'Genotyped', 'ALT_Frq', 'MAF', 'Rsq']):
    # Only need these columns to merge
    i = 0
    df_merged = '' # A DataFrame of merged .info.gz df
    lst_index_col_names = [] # A list to store column names of index columns

    # Create POS column to merge on, since SNP in different files may have flipped REF and ALT, but POS is always the same
    # Convert MAF, Alt_frq and Rsq to numeric
    col_name_alt_frq, col_name_maf, col_name_r2 = cols_to_keep[-3:]
    for df in lst_info_df:
        df['POS'] = df['SNP'].apply(lambda x: int(x.split(':')[1]))
        df[col_name_alt_frq] = pd.to_numeric(df[col_name_alt_frq], errors='coerce')
        df[col_name_maf] = pd.to_numeric(df[col_name_maf], errors='coerce')
        df[col_name_r2] = pd.to_numeric(df[col_name_r2], errors='coerce')

    while i<len(lst_info_df):
        if i==0: # Merge the 1st and 2nd df
            # Need to merge on multiple keys: ['SNP', 'REF(0)', 'ALT(1)','Genotyped']
            # Otherwise will be NAs if a variant is not found in the first info file
            # Use reset_index() to track index number of a variant in each input file
            # Eg. variant rs0000 is from row index 1 in input file 1 and row #3 in input file 2
            df_merged = lst_info_df[i][cols_to_keep+['POS']].reset_index().merge(lst_info_df[i+1][cols_to_keep+['POS']].reset_index(),
                                                           how='outer',
                                                           on=['SNP', 'POS'],
                                                           suffixes=('_group'+str(i+1), '_group'+str(i+2)))

            lst_index_col_names.append('index_group'+str(i+1))
            lst_index_col_names.append('index_group'+str(i+2))
            i = i + 2
        else:
            # Merge and rename columns of df with correct suffix (ie. _groupX)
            df_merged = df_merged.merge(lst_info_df[i][cols_to_keep+['POS']].reset_index().rename(columns={'MAF':'MAF_group'+str(i+1),
                                                                                                           'Rsq':'Rsq_group'+str(i+1),
                                                                                                           'ALT_Frq':'ALT_Frq_group'+str(i+1),
                                                                                                           'Genotyped':'Genotyped_group' + str(i+1),
                                                                                                           'index':'index_group' + str(i+1)}),
                                        how='outer',
                                        on=['SNP', 'POS'])
            lst_index_col_names.append('index_group' + str(i+1))
            i += 1
    # Get chr from snp IDs, maybe add this feature in the future
    # df_merged['chr'] = df_merged['SNP'].apply(lambda x: x.split(':')[0])
    return df_merged, lst_index_col_names


# This function determines genotype status of each variant
# Patameters:
# - df: a DataFrame to be processes (df_merge in this script)
# - genotyped_label: label(s) for genotyped variants ('TYPED' and 'TYPED_ONLY' in TOPmed)
# - imputed_label: label(s) for imputed variants ('IMPUTED' in TOPmed)
# - number_of_input_files: Number of input files
# Return:
# - df will be modified inplace such that
#   - Genotype status of individual input files (Genotyped_group1, Genotyped_group2, .. etc.) will be removed
#   - A column named 'Genotyped' of final mixed genotype status will be retained
#       - If mixed_genotype=True, final genotype status are: ALL=All genotyped; SOME=Some genotyped; NONE=None genotyped
#       - If mixed_genotype=False, final genotype status will be the same as in first input file, or the first file a given SNP is found.
#       - If some *.info.gz files contain ALL/SOME/NONE while others use GENOTYPED/TYPED, final genotype status will be mixed of ALL/SOME/NONE/GENOTYPED/TYPED
def get_genotype_status(df, mixed_genotype, genotyped_labels, imputed_labels, number_of_input_files):
    # Get column names of genotype status of each input file
    lst_genotype_col_names = [f'Genotyped_group{i+1}' for i in range(number_of_input_files)]
    if not mixed_genotype: # Use the first genotype status available
        df['Genotyped'] = df[lst_genotype_col_names].fillna(method='backfill', axis=1)[lst_genotype_col_names[0]]
    else:
        lst_genotyped_label = genotyped_labels.split('/')
        lst_imputed_label = imputed_labels.split('/')
        print('\nGenotype status label(s) used for genotyped variants:', lst_genotyped_label)
        print('Genotype status label(s) used for imputed variants:', lst_imputed_label)

        df['Genotyped'] = 'NONE' # Prefill with 'NONE'
        # Find ALL and SOME variants, change their genotype status
        mask_total_gped = df[lst_genotype_col_names]==lst_genotyped_label[0]
        for gp_label in lst_genotyped_label[1:]+['ALL']: # Include 'ALL' as label to 'genotyped' variants
            mask_total_gped = mask_total_gped | (df[lst_genotype_col_names]==gp_label)
        mask_all = mask_total_gped.all(axis=1)

        # If genotype status of any input file is 'SOME', a given variant is also 'SOME'
        # Find variants with at least one 'genotyped' status and fill them with 'SOME', then reassign SNPs with only genotyped status with 'ALL'
        mask_any_gped = (mask_total_gped.any(axis=1)) | (df[lst_genotype_col_names]=='SOME').any(axis=1)

        df.loc[mask_any_gped, 'Genotyped'] = 'SOME'
        df.loc[mask_all, 'Genotyped'] = 'ALL'
    df.drop(columns=lst_genotype_col_names, inplace=True)


# This function calculated r2 (use method defined by dict_flags['--r2_output'])， MAF, and altFrq,
# assign values to a column (col_name_r2_combined, col_name_maf_combined) in df_merged
# Parameters:
# - df_merged: merged dataframe from all .info.gz files
# - col_name_r2_combined: column name of combined r2
# - col_name_maf_combined: column name of combined minor allele frequency
# - col_name_alt_frq: column name of combined alternative allele frequency (ALT_Frq)
# - dict_flags: argument dictionary of all flags
# (Deprecated) - lst_input_fn: a list of names of input files
# (Deprecated) - r2_merge_type: use value of dict_flags['--r2_output'], can be 'weighted_average', 'first', 'mean', 'min' and 'max'
#                  Defines how Rsq is calculated
def __calculate_r2_maf_altFrq(df_merged, col_name_r2_combined,
                              col_name_maf_combined, col_name_alt_frq_combined, dict_flags):
    lst_number_of_individuals = []  # A list to store number of individuals of each input file
    # Count number of individuals in each input file
    for fn in dict_flags['--input']:
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
    lst_col_r2_names=[] # Store column names of r2 from each file, such as Rsq_group1, Rsq_group2 ...
    for i in range(len(dict_flags['--input'])):
        col_name_r2 = 'Rsq_group'+str(i+1)
        col_name_maf = 'MAF_group' + str(i+1)
        col_name_alt_frq = 'ALT_Frq_group'+str(i+1)
        lst_col_r2_names.append(col_name_r2)

        # !!! Must use pd.DataFrame.sum(), mean(), etc. functions to deal with missing values (NaN)
        # Otherwise will get NaN if adding NaN to another value directly.
        df_merged[col_name_r2+'_adj'] = df_merged[col_name_r2] * lst_number_of_individuals[i] # Rsq * # individuals
        df_merged[col_name_maf+'_adj'] = df_merged[col_name_maf] * lst_number_of_individuals[i]   # MAF * # individuals
        df_merged[col_name_alt_frq + '_adj'] = df_merged[col_name_alt_frq] * lst_number_of_individuals[i]  # MAF * # individuals

        # Weight is number of individuals (but not include missing values)
        # Code below reflects missing values by * and / the same values (Rsq)
        # Missing values are not counted when calculate weighted Rsq, MAF and AF (so use if for all!)
        df_merged[col_name_r2 + '_weight_adj'] = lst_number_of_individuals[i] * df_merged[col_name_r2]/ df_merged[col_name_r2]
        lst_col_names_r2_adj.append(col_name_r2+'_adj')
        lst_col_names_maf_adj.append(col_name_maf + '_adj')
        lst_col_names_alt_frq_adj.append(col_name_alt_frq + '_adj')
        lst_col_names_weight_adj.append(col_name_r2 + '_weight_adj')

    df_merged[col_name_alt_frq_combined] = df_merged[lst_col_names_alt_frq_adj].sum(axis=1) / df_merged[lst_col_names_weight_adj].sum(axis=1)
    df_merged[col_name_maf_combined] = df_merged[lst_col_names_maf_adj].sum(axis=1) / df_merged[lst_col_names_weight_adj].sum(axis=1)
    if dict_flags['--r2_output'] == 'weighted_average': # ie. dict_flags['--r2_output']=='weighted_average'
        df_merged[col_name_r2_combined] = df_merged[lst_col_names_r2_adj].sum(axis=1) / df_merged[lst_col_names_weight_adj].sum(axis=1)
    elif dict_flags['--r2_output'] == 'z_transformation': # Fisher's z-transformation
        # Fisher's z-transformation -> weighted average -> tanh
        # z transformation: 0.5 * np.log((1 + r2) / (1 - r2))
        lst_col_names_r2_z_trans = []
        lst_col_names_r2_z_trans_weight_adj = []
        for i in range(len(dict_flags['--input'])):
            col_name_r2 = 'Rsq_group'+str(i+1)
            # adjust with --r2_offset when Rsq=1 (Rsq-1). This column will be saved in output
            df_merged['Rsq_group'+str(i+1)+'_offset_adj']=df_merged['Rsq_group'+str(i+1)]
            mask_rsq = (df_merged['Rsq_group'+str(i+1)+'_offset_adj']==1)
            # Substract --r2_offset when R2=1
            df_merged.loc[mask_rsq, 'Rsq_group' + str(i + 1) + '_offset_adj'] -= dict_flags['--r2_offset']

            # Cannot use apply(), it will raise 'division by zero error' when Rsq=1
            # df_merged[col_name_r2+'_z_trans'] = df_merged['Rsq_group'+str(i+1)].apply(lambda x: 0.5 * np.log((1 + x) / (1 - x)))
            # The code below will return np.inf when Rsq=1. Then np.tanh(np.inf)=1
            df_merged[col_name_r2 + '_z_trans'] = 0.5 * np.log((1 + np.sqrt(df_merged['Rsq_group'+str(i + 1)+'_offset_adj']) ) / (1 - np.sqrt(df_merged['Rsq_group'+str(i + 1)+'_offset_adj']) ))
            df_merged[col_name_r2 + '_z_trans_weight_adj'] = df_merged[col_name_r2+'_z_trans'] * lst_number_of_individuals[i] # Z_transed_r * weight
            lst_col_names_r2_z_trans.append(col_name_r2 + '_z_trans')
            lst_col_names_r2_z_trans_weight_adj.append(col_name_r2 + '_z_trans_weight_adj')
        df_merged['R_z_trans_combined'] = df_merged[lst_col_names_r2_z_trans_weight_adj].sum(axis=1) / df_merged[lst_col_names_weight_adj].sum(axis=1)
        df_merged[col_name_r2_combined] = (np.tanh(df_merged['R_z_trans_combined']))**2
        df_merged['var_of_z_trans_r'] = df_merged[lst_col_names_r2_z_trans].var(axis=1) # variance of z transformed R
    elif dict_flags['--r2_output'] == 'first': # ie. dict_flags['--r2_output']=='first'
        df_merged[col_name_r2_combined] = df_merged['Rsq_group1']
    elif dict_flags['--r2_output'] == 'mean':
        df_merged[col_name_r2_combined] = df_merged[lst_col_r2_names].mean(axis=1)
    elif dict_flags['--r2_output'] == 'min':  # ie. dict_flags['--r2_output']=='min'
        df_merged[col_name_r2_combined] = df_merged[lst_col_r2_names].min(axis=1)
    else:  # ie. dict_flags['--r2_output']=='max'
        df_merged[col_name_r2_combined] = df_merged[lst_col_r2_names].max(axis=1)

    # Drop columns that are not needed in output
    df_merged.drop(labels=lst_col_names_r2_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_maf_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_alt_frq_adj, axis=1, inplace=True)
    df_merged.drop(labels=lst_col_names_weight_adj, axis=1, inplace=True)
    return lst_number_of_individuals    # Return number of individuals in each file


# This function save excluded and kept variants into .txt files
# Parameters:
#  - df_merged: A DataFrame of variants from all input files with desired fields
#  - missing: A variant is allowed to be missing in this number of files. Otherwise it will be excluded
#  - lst_index_col_names: original indices of each input dataframe are stored in these column for sorting purpose
# Output:
#  - Two text files saved in current directory: variants_kept.txt and variants_excluded.txt
def __process_output(df_merged, dict_flags, lst_index_col_names):
    # Set output file names
    to_keep_fn = dict_flags['--output']+'_variants_retained.info.txt'
    to_exclude_fn = dict_flags['--output']+'_variants_excluded.info.txt'
    mask_to_keep, mask_to_exclude = '', ''   # Use this when --missing is 0

    lst_col_names = []  # A list to store column names of Rsq (such as Rsq_group1, Rsq_group2, etc.)
    lst_genotype_status = []
    for i in range(len(dict_flags['--input'])):
        col_name_r2 = 'Rsq_group' + str(i + 1)
        # col_name_genotype = 'Genotyped_group' + str(i + 1) # Genotype status
        # lst_genotype_status.append(col_name_genotype)

        # col_name_maf = 'MAF_group' + str(i + 1)
        # col_name_alt_frq = 'ALT_Frq_group' + str(i + 1) # Column name of alternative allele
        # # Convert values of r2, MAF and ALT_Frq columns from strings to numbers for later calculations
        # df_merged[col_name_r2] = pd.to_numeric(df_merged[col_name_r2], errors='coerce')
        # df_merged[col_name_maf] = pd.to_numeric(df_merged[col_name_maf], errors='coerce')
        # df_merged[col_name_alt_frq] = pd.to_numeric(df_merged[col_name_alt_frq], errors='coerce')

        # Process retained and excluded SNPs
        if dict_flags['--missing'] == 0:  # If only save variants shared by all input files
            # Check columns of Rsq (MAF also works) such as: Rsq_group1, Rsq_group2,...
            # If MAF_groupX is NaN, then this variant is missing from the corresponding groupX.
            if i == 0:  # If the first iteration
                mask_to_keep = df_merged[col_name_r2].notna()
                mask_to_exclude = df_merged[col_name_r2].isna()
            else: # Add additional mask in each iteration
                mask_to_keep = mask_to_keep & df_merged[col_name_r2].notna()
                mask_to_exclude = mask_to_exclude | df_merged[col_name_r2].isna()

            # Check if need to filter for SNPs using given r2 threshold
            if dict_flags['--r2_threshold'] != 0:
                mask_to_keep = mask_to_keep & (df_merged[col_name_r2] >= dict_flags['--r2_threshold'])
                mask_to_exclude = mask_to_exclude | (df_merged[col_name_r2] < dict_flags['--r2_threshold'])

        lst_col_names.append(col_name_r2)  # Save column names of Rsq_group1, Rsq_group2, etc.

    # Determine genotype status: ALL=All genotyped, SOME=Some genotype, NONE=none genotyped
    get_genotype_status(df_merged, dict_flags['--mixed_genotype_status'], dict_flags['--genotyped_label'],
                        dict_flags['--imputed_label'], len(dict_flags['--input']))

    # Calculate combined r2 (first, weighted_average or mean)
    col_name_r2_combined = 'Rsq_combined'
    col_name_maf_combined = 'MAF_combined'
    col_name_alt_frq_combined = 'ALT_Frq_combined'

    lst_number_of_individuals = __calculate_r2_maf_altFrq(df_merged, col_name_r2_combined,
                                                          col_name_maf_combined, col_name_alt_frq_combined, dict_flags)

    # Save result to output files
    float_format = '%.6f'
    if dict_flags['--missing'] == 0:
        # Use %g for precision formatting, will get a mix of scientific notation when number is too small
        # Or use %.6f, for 6 precision after dot. But if value is stoo small will be round to 0
        # Save to *.variant_retained.info.txt sorted by position and index in each input file,
        # since there could be multiple variants at the same position
        df_merged[mask_to_keep].sort_values(by=['POS']+lst_index_col_names)\
            .drop(columns=lst_index_col_names+['POS'])\
            .to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'], float_format=float_format)
        df_merged[mask_to_exclude].drop(columns=lst_index_col_names+['POS']).to_csv(to_exclude_fn, index=False,
                                                                                           sep='\t', na_rep=dict_flags['--na_rep'],
                                                                                           float_format=float_format)
        print('\tNumber of saved variants:', len(df_merged[mask_to_keep]))
        print('\tNumber of excluded variants:', len(df_merged[mask_to_exclude]))

        # Write important info into log file
        log_fn = dict_flags['--output'] + '.log'  # Save Important processing info into a .log file for user reference
        with open(log_fn, 'a') as log_fh:
            log_fh.write('\tTotal number of all input files combined: '+str(len(df_merged))+'\n')
            log_fh.write('\tNumber of saved variants:'+str(len(df_merged[mask_to_keep]))+'\n')
            log_fh.write('\tNumber of excluded variants:'+str(len(df_merged[mask_to_exclude]))+'\n')
    else:   # Use user specify allowed missingness, max value is number of input files
        # Calculate number of missing files (by number of NA value of Rsq column)
        df_merged['missing'] = df_merged[lst_col_names].isna().sum(axis=1)

        # Filter by dict_flags['--r2_threshold'], if merged Rsq>=dict_flags['--r2_threshold'] then keep this variant
        # Also filter by number of missing files (dict_flags['--missing'])
        mask_to_keep = (df_merged['missing']<=dict_flags['--missing']) & \
                       (df_merged[col_name_r2_combined]>=dict_flags['--r2_threshold'])
        mask_to_exclude = (df_merged['missing']>dict_flags['--missing']) | \
                          (df_merged[col_name_r2_combined]<dict_flags['--r2_threshold'])

        df_merged.loc[mask_to_keep].sort_values(by=['POS']+lst_index_col_names)\
                .drop(columns=lst_index_col_names+['POS', 'missing'])\
                .to_csv(to_keep_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'], float_format=float_format)
        df_merged.loc[mask_to_exclude].drop(columns=lst_index_col_names+['POS', 'missing'])\
            .to_csv(to_exclude_fn, index=False, sep='\t', na_rep=dict_flags['--na_rep'], float_format=float_format)

        print('\tNumber of saved variants:', len(df_merged.loc[mask_to_keep]))
        print('\tNumber of excluded variants:', len(df_merged.loc[mask_to_exclude]))

        # Write important info into log file
        log_fn = dict_flags['--output'] + '.log'  # Save Important processing info into a .log file for user reference
        with open(log_fn, 'a') as log_fh:
            log_fh.write('\tTotal number of all input files combined: '+str(len(df_merged))+'\n')
            log_fh.write('\tNumber of saved variants: '+str(len(df_merged.loc[mask_to_keep]))+'\n')
            log_fh.write('\tNumber of excluded variants: '+str(len(df_merged.loc[mask_to_exclude]))+'\n')

    # *.index.text is for debugging purpose, it is not needed in merge
    if dict_flags['--verbose']:
        df_merged.sort_values(by=['POS']+lst_index_col_names)[['SNP']+lst_index_col_names].to_csv(dict_flags['--output']+'_index.txt',
                                                                                              sep='\t', index=False, na_rep='-')
    print('\nNumbers of individuals in each input file:', lst_number_of_individuals)

    # Write into a log file
    log_fn = dict_flags['--output'] + '.log'  # Save Important processing info into a .log file for user reference
    with open(log_fn, 'a') as log_fh:
        log_fh.write('\nNumbers of individuals in each input file: '+str(lst_number_of_individuals))

    return lst_number_of_individuals, len(df_merged[mask_to_keep]) # Return list of number of individuals and number of SNPs kept
# ---------------- End opf helper functions -----------------

# A wrapper function to run this script
# Returns dict_flags for next step
def get_snp_list(dict_flags):
    lst_info_df = __get_lst_info_df(dict_flags)
    df_merged, lst_index_col_names = __merge_snps(lst_info_df)
    print('\nNumber of variants:')
    print('\tTotal number from all input files:', len(df_merged))
    lst_number_of_individuals, number_snps_kept = __process_output(df_merged, dict_flags, lst_index_col_names)
    return lst_number_of_individuals, number_snps_kept # These return values are used in the final merging step

