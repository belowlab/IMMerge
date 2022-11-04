# Step 02
# This code is used to merge imputed files from TopMed based on list of saved variants in prior steps
# The goal is to combine post-imputation vcf.gz files into one vcf.gz file

from xopen import PipedCompressionWriter, xopen # Faster than gzip
import os
# import process_args
# from .process_args import process_args
# from .get_SNP_list import get_snp_list
# from .check_r2_setting_for_imputation import check_imputatilson_parameters
import time

if __name__ == "__main__":
    from process_args import process_args
    from get_SNP_list import get_snp_list
else: # import this way when used as a module. '.' indicates to import from the current module (directory)
    from .process_args import process_args
    from .get_SNP_list import get_snp_list

# ################################ Helper functions ################################
# Print progress par in console
# - progress: current progress (number of SNPs processed)
def progress_bar(progres):
    total = number_snps_kept # total number of SNPs needs to be processed
    percent = 100 * (progres/total)
    bar = '=' * int(percent) + '-' * int(100 - percent)
    print(f'|{bar}| {percent:.2f}%', end='\r')


# This function prints execution time at the end of run
def print_execution_time(start_time):
    duration = time.time() - start_time
    hours = duration // 3600
    duration = duration % 3600
    minutes = duration // 60
    seconds = duration % 60
    print('\nDuration: {:.2f} hours, {:.2f} minutes, {:.2f} seconds'.format(hours, minutes, seconds))
    print('End of run')

    global LOG_TXT
    LOG_TXT += f'\n\nDuration: {hours:.2f} hours, {minutes:.2f} minutes, {seconds:.2f} seconds'
    LOG_TXT += '\n\nEnd of run'


# This function writes IMMerge arguments in to meta information of merged output
# Lines start with "##"
def write_args(fh_output):
    comment_line = '##IMMerge_command=merge_files.py '
    for k, v in dict_flags.items():
        if k=='--input' or k=='--info':
            comment_line += f'{k}'
            for fn in v:
                comment_line += ' ' + fn
        else:
            comment_line += f' {k} {v}'
    comment_line = comment_line + '; Date=' + time.strftime('%a %b %d %H:%M:%S %Y', time.localtime()) + '\n'
    fh_output.write(comment_line)

# Check if there are duplicated samples, then rename duplicate samples as ID, ID:2, ID:3,...
# Output duplicate IDs to file dup_IDs.txt
# Return merged column header to write into output file
def rename_duplicated_samples(lst_column_header, inx_indiv_id_starts):
    path, _ = os.path.split(dict_flags['--output']) # Get output directory
    dup_id_fn = 'dup_IDs.txt' # Write duplicate ID s to dup_IDs.txt in output directory
    dup_id_fh = open(os.path.join(path, dup_id_fn), 'w')
    print('\tChecking sample ID duplication: Duplicate sample IDs are saved in', os.path.join(path, dup_id_fn))
    global LOG_TXT
    LOG_TXT += '\tChecking sample ID duplication: Duplicate sample IDs are saved in '+os.path.join(path, dup_id_fn) + '\n'

    merged_header = lst_column_header[0] # merged column header to return
    # dict_id_count is a dictionary to track number after colon to rename an ID
    dict_id_count = {k:1 for k in lst_column_header[0].split(maxsplit=inx_indiv_id_starts)[-1].split()}
    for i in range(1, len(lst_column_header)): # Check headers from non-first input files
        lst_tmp = lst_column_header[i].split(maxsplit=dict_flags['--duplicate_id']+inx_indiv_id_starts)[-1].split()
        for j in range(len(lst_tmp)):
            if dict_id_count.get(lst_tmp[j]) is None: # If current ID is not seen before
                dict_id_count[lst_tmp[j]] = 1
            else:
                dict_id_count[lst_tmp[j]] += 1
                dup_id_fh.write(lst_tmp[j] + ':' + str(dict_id_count[lst_tmp[j]])+'\n')
                lst_tmp[j] = lst_tmp[j] + ':' + str(dict_id_count[lst_tmp[j]]) # Replace current ID with ID:index+1

        merged_header += '\t' + '\t'.join(lst_tmp) # Reconstruct the header line
    dup_id_fh.close()
    return merged_header

# Merge header lines of input files (.dose.vcf.gz) together and write into a BGZF file (bgzip file)
# Parameters:
# - fh_output_raw: file handle of output file (need to use bgzip module for actual writing)
# - lst_input_fh: list of file handles of input files
def merge_header_lines(lst_input_fh, fh_output):
    global LOG_TXT
    # Read trough info lines (line start with ##) and write into output file
    # Process all files at the same time. DO NOT assume All input files have the same number of meta info lines
    lst_column_header = [] # Store lines for column headers of each input files
    if dict_flags['--meta_info'] == 'none': # Do not include any meta information
        for fh in lst_input_fh:
            # Read lines until reach column header line (starts with "#")
            line = fh.readline().strip()
            while line[:2] == '##':
                line = fh.readline().strip()
            lst_column_header.append(line)
    elif dict_flags['--meta_info'] == 'all': # Include all meta information
        for fh in lst_input_fh:
            # Read lines until reach column header line (starts with "#"), and write to output
            line = fh.readline().strip()
            while line[:2] == '##':
                fh_output.write(line+'\n')
                line = fh.readline().strip()
            lst_column_header.append(line)
    else: # Only include meta information from a given file by index (1-based)
        for i in range(len(lst_input_fh)):
            if i+1 == dict_flags['--meta_info']:
                line = lst_input_fh[i].readline().strip()
                while line[:2] == '##':
                    fh_output.write(line + '\n')
                    line = lst_input_fh[i].readline().strip()
            else: # Ignore meta info from other input files
                line = lst_input_fh[i].readline().strip()
                while line[:2] == '##':
                    line = lst_input_fh[i].readline().strip()
            lst_column_header.append(line)

    # Write IMMerge arguments after handling meta information
    write_args(fh_output)

    # Merge column headers
    if line[0] != '#': # Sanity check, can also use: assert line[0] == '#', 'header line should start with #'
        print('Error: header line should start with #\nExit')
        LOG_TXT += 'Error: header line should start with #\nExit\n'
        exit()

    # In current version of VCF format (4.3) these are fixed, mandatory columns, then followed by individual IDs:
    #   - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    inx_indiv_id_starts = line.split(maxsplit=10).index('FORMAT') + 1  # Find from which column genotype data starts
    inx_info_column = line.split(maxsplit=10).index('INFO') # Find index of column "INFO"
    print('\tIndividual genotype data starts from column #:', inx_indiv_id_starts + 1)
    LOG_TXT += f'\tIndividual genotype data starts from column #: {inx_indiv_id_starts + 1}\n'

    if dict_flags['--check_duplicate_id']: # Check if there are duplicated samples and rename
        line_to_write = rename_duplicated_samples(lst_column_header, inx_indiv_id_starts)
    else:
        # Create merged header line, initialize with header from the first file
        line_to_write = lst_column_header[0]
        for col_header in lst_column_header[1:]:  # Process header line of remaining input files
            # Set maxsplit in slit() function to remove duplicated IDs
            line_to_write = line_to_write + '\t' + \
                            col_header.split(maxsplit=dict_flags['--duplicate_id'] + inx_indiv_id_starts)[-1]
    fh_output.write((line_to_write + '\n'))
    return inx_indiv_id_starts, inx_info_column


# This function read line of a given file handle, until find line of given SNP (by parameter snp))
# Return string of that line (if not found then return an empty string '')
def search_SNP_and_read_lines(snp, fh):
    inx_snp_ID = 2  # Index of column that contains variant ID such as chr21:10000777:C:A (index is hard coded as this column is fixed in VCF format)
    line = fh.readline().strip()  # Read in each line to search for the given variant
    # input_snp = line.split(maxsplit=inx_snp_ID + 1)[-2]
    while line != '':  # Keep read line until find given snp
        if dict_flags['--use_rsid']: # If ID column contains rsID instead of chr:pos:ref:alt, need to construct the right format for comparisons
            tmp_lst = line.split(maxsplit=10) # First 8 columns are fixed, so it is safe to set maxsplit to 10
            # Indices of #CHROM, POS, REF and ALT are 0, 1, 3 and 4
            input_snp = f'{tmp_lst[0]}:{tmp_lst[1]}:{tmp_lst[3]}:{tmp_lst[4]}'
        else:
            input_snp = line.split(maxsplit=inx_snp_ID + 1)[inx_snp_ID]

        if input_snp == snp: # if found
            break
        else:
            line = fh.readline().strip()  # Read in each line to search for the given variant

    if line == '': # If reach end of file but the given SNP is not found, something is wrong
        print(f'\nError: reached end of file {fh.name}, but SNP {snp} is not found.')
        print('- Do your *.info.gz files contain extra SNPs not found in VCFs? Try using make_info.py to create new info files\nExit')
        exit()
    return line

# This function takes in a variant ID,
# checks every input file to find corresponding line and returns a merged line.
# If a snp is not found in certain input file, fill with symbol of missing values for that file
# Missing values symbol can be any user supplied string, default is '.'
# Value of INFO field is replace with new values calculated from all input files
# Parameters:
# - snp:
# - number_of_dup_id: dict_flag['--duplicate_id']
# - inx_info_column: index of INFO columns in input .dose.vcf.gz file
# - new_info_val: new info values to replace value of the INFO column in input .dose.vcf.gz file
# - inx_indiv_id_starts: index of the first individual ID (such as R200000348) columns in input .dose.vcf.gz file
# - lst_alt_frq_val: a list of alt_Frq values from all input files. Use this to decide if a SNP is missing from a input file
# - lst_input_fh: list of input file handles
def merge_individual_variant(snp, number_of_dup_id, inx_info_column, new_info_val, inx_indiv_id_starts, lst_alt_frq_val, lst_input_fh):
    merged_line = ''
    # flag_first_na: An indicator to show whether a variant is missing from the first input file
    # If missing, need to fill info columns from other input file
    flag_first_na = False
    for i in range(len(lst_input_fh)):
        if i == 0:  # If this is the first input file
            if lst_alt_frq_val[i] != dict_flags['--na_rep']:
                # If this is the first input file and variant is not missing,
                # then keep all columns including info columns
                lst_tmp = search_SNP_and_read_lines(snp, lst_input_fh[i]).split(maxsplit=inx_info_column+1)

                lst_tmp[inx_info_column] = new_info_val
                merged_line = '\t'.join(lst_tmp)
            else:
                # If variant is missing from the file, then fill with '.|.' (missing genotype) and exclude the first few info columns
                # Need to get values of info columns from other input files
                # merged_line = '\t'.join([dict_flags['--na_rep']+'|'+dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i])])
                merged_line = '\t'.join([dict_flags['--na_rep'] + '|' + dict_flags['--na_rep'] for j in
                                         range(lst_number_of_individuals[i])])
                flag_first_na = True  # Indicate flag: Fill info columns later with other input files
        else:  # If this is not the first file
            if lst_alt_frq_val[i] == dict_flags['--na_rep']:
                # If variant is missing from this input file, then fill 'NA' or user provided missing value symbol
                merged_line = merged_line + '\t' + '\t'.join(
                    [dict_flags['--na_rep']+'|'+dict_flags['--na_rep'] for j in range(lst_number_of_individuals[i]-number_of_dup_id)])
            else:  # If not missing, split and remove info columns
                line = search_SNP_and_read_lines(snp, lst_input_fh[i])
                # Skip number_of_dup_id samples when merging
                merged_line = merged_line + '\t' + line.split(maxsplit=number_of_dup_id+inx_indiv_id_starts)[-1]
                if flag_first_na:  # If info columns is still missing, add these columns back
                    lst_tmp = line.split(maxsplit=inx_indiv_id_starts)[:-1]
                    lst_tmp[inx_info_column] = new_info_val # replace INFO field with updated INFO value
                    info_cols = '\t'.join(lst_tmp)
                    merged_line = info_cols + '\t' + merged_line
                    flag_first_na = False
    return merged_line

# This function reads in variant_kept.txt file as a reference, and merge rest of lines in input files
# (Header lines should be handled by merge_header_lines() function, so they are ignored in this function)
# Merged lines are written into output file
def merge_files(dict_flags, inx_info_column, inx_indiv_id_starts, lst_input_fh, fh_output):
    global LOG_TXT
    # Load variants_kept.txt as a reference of merged file
    # Then compare snp ids in each input file with the SNPs in variants_kept.txt
    # And replace ALT_frq, MAF and r2 values in INFO column with new values calculated from all input files
    # - ie, use these columns in variant_kept.txt file: ALT_Frq_combined, MAF_combined, Rsq_combined, Genotyped
    count = 0  # For console output
    fh_snp_kept = open(dict_flags['--output']+'_variants_retained.info.txt')
    print('\t'+dict_flags['--output']+'_variants_retained.info.txt file loaded\n')
    LOG_TXT += '\t'+dict_flags['--output']+'_variants_retained.info.txt file loaded\n\n'

    # Read the first line (header), find index of combined ALT_freq, MAF, Rsq
    line = fh_snp_kept.readline().strip().split()
    indx_ALT_Frq_combined = line.index('ALT_Frq_combined')
    indx_MAF_combined = line.index('MAF_combined')
    indx_Rsq_combined = line.index('Rsq_combined')
    indx_genotyped = line.index('Genotyped')

    while True:
        line_snp_kept = fh_snp_kept.readline().strip()
        if line_snp_kept != '':
            lst_snp_kept = line_snp_kept.split()
            snp = lst_snp_kept[0]  # Get ID of SNP to be kept

            # Get ALT_frq, MAF and r2 values
            # Use ALT_frq_groupX values to decide whether SNP is missing in the Xth file
            lst_alt_frq_val = []
            ALT_Frq_combined = lst_snp_kept[indx_ALT_Frq_combined]
            MAF_combined = lst_snp_kept[indx_MAF_combined]
            Rsq_combined = lst_snp_kept[indx_Rsq_combined]
            genotyped = lst_snp_kept[indx_genotyped]
            # Create new INFO value to replace INFO column in merged line
            # new_info_val = 'AF='+ALT_Frq_combined+';MAF='+MAF_combined+';R2='+Rsq_combined+';'+genotyped
            new_info_val = f'AF={ALT_Frq_combined};MAF={MAF_combined};R2={Rsq_combined};{genotyped}'

            inx_alt_frq_col = 4  # Index of ALT_frq column group1 (since variant_kept file is build by the program, the index can be hard coded)
            # Add ALT_frq of group2,...,groupX to each list
            for i in range(len(lst_input_fh)):
                lst_alt_frq_val.append(lst_snp_kept[inx_alt_frq_col])
                inx_alt_frq_col += 3

            # Search line of this snp in each input file, fill with missing values if not exist
            # Use ALT_Frq_groupX column in variant_kept.txt file to indicate if this snp exist in input files
            # Eg. if ALT_Frq_group1 is missing, then this snp does not exist in the first input file
            # print('SNP to be kept:', snp)
            merged_line = merge_individual_variant(snp, dict_flags['--duplicate_id'], inx_info_column, new_info_val,
                                                   inx_indiv_id_starts, lst_alt_frq_val, lst_input_fh)
            # Write merged line to output file
            fh_output.write(merged_line)
            fh_output.write('\n')

            # Keep console busy
            count += 1
            if count%100 == 0: progress_bar(count)
        else:
            bar = '=' * 100 # Complete progress bar, fill with 100%
            print(f'|{bar}| 100.00%', end='\r')
            print('\nEnd of snp to be kept file, total of {} variants merged'.format(count))
            break  # If reach end of SNP to be kept, then stop
    # Close file handles when done
    fh_snp_kept.close()


# ################################ End of helper functions ################################



def run_merge_files(args):
    '''
    A function to run merging steps
    Params:
    - args: a dictionary containing arguments when IMMerged is used as a python module, formatted as {'--input':['file1', 'file2'], '--missing':1}. Valid flags are:
        --input: (Required) Files to be merged, multiple files are allowed. Must in gzipped or bgziped VCF format.
        --info: (Optional) Directory/name to info files.
                Default path is the same directory as corresponding input file, default info files share the same name as input file, except for suffix (.info.gz).
        --output: (Optional) Default is "merged.vcf.gz" and saved at current working directory. Output file name without suffix.
        --thread: (Optional) Default value is 1. Defines how many thread to use in multiprocessing. If number of threads <0, will use 1 instead of user supplied value.
        --missing: (Optional) Default is 0. Defines number of missing values allowed for each variant.
        --na_rep: (Optional) Default is "." (ie. ".|." for genotype values). Defines what symbol to use for missing values. This flag is ignored if --missing is 0.
        --r2_threshold: (Optional) Default is 0, ie. no filtering. Only variants with combined imputation quality score r2≥r2_threshold will be saved in the merged file.
        --r2_output: (Optional) Default is "z_transformation". Defines how imputation quality score is calculated in the output file.
        --r2_offset: (Optional) Default is 0.001. Adjust R squared by --r2_offset if imputation quality Rsq=1. Only valid for z transformation to avoid infinity.
        --duplicate_id: (Optional) Default is 0. Defines number of duplicate individuals in each input file. Duplicated IDs should be the first N columns in each file.
        --check_duplicate_id: (Optional) Default is False. Check if there are duplicate IDs, then rename non-first IDs to ID:2, ID:3, ..., ID:index_of_input_file+1.
        --write_with: (Optional) Default is bgzip. Write to bgziped file with bgzip. User can supply specific path to bgzip such as /user/bin/bgzip.
        --meta_info: (Optional) Valid values are {index of input file (1-based), 'none', 'all'}.
                    What meta information (lines start with '##') to include in output file. Default is 1 (meta information from the first input file).
        --use_rsid: (Optional) Default is False. If input VCFs use rsID instead of chr:pos:ref:alt, set this option to True to avoid duplicate IDs (rsID may not be unique).
                    New IDs in chr:pos:ref:alt format instead of rsID will be used to merge.
                    Use make_info.py to make info files, or follow the required format.
        --verbose: (Optional) Default is False. Print more messages.
    Return:
    - Save a merged and bgziped vcf file as specified in --output
    '''
    start_time = time.time()  # Track execution time

    global VERSION  # Global variable Version of IMMerge
    VERSION = '0.0.4'
    global LOG_TXT
    LOG_TXT = ''  # Global variable to track info printed. Write to log file when run is finished

    print(f'\nIMMerge version {VERSION}')
    print('Job starts at (local time) ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))
    LOG_TXT += f'IMMerge Version {VERSION}\n'
    LOG_TXT += 'Job starts at (local time) ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)) + '\n'

    # Process terminal input or arguments passed through python module call
    # dict_flags contains values for IMMerge flags:
    global dict_flags

    if args != '': # If used as a module
        # Parse arguments, convert to a list and pass to process_args()
        args_list = []
        for k, v in args.items():
            args_list.append(k)
            if isinstance(v, list): # add all values if it is a list
                for val in v: args_list.append(val)
            else: args_list.append(str(v)) # Args from command line are strings, so need to convert for process_args() to works
        dict_flags = process_args(args_list)
    else: # If called from terminal, so args==''
        dict_flags = process_args(args)  # Process terminal input, use it as a global variable

    for k, v in dict_flags.items():
        LOG_TXT += '\t' + k + ': ' + str(v) + '\n'

    # Do not need this one below
    # check_r2_setting_for_imputation.check_imputation_parameters(lst_fn=dict_flags['--input'])

    global lst_number_of_individuals
    global number_snps_kept
    lst_number_of_individuals, number_snps_kept = get_snp_list(dict_flags)
    LOG_TXT += f'\nNumber of variants retained: {number_snps_kept}\n'
    LOG_TXT += f'Numbers of individuals in each input file: {lst_number_of_individuals}\n'

    print('\nStart merging files:')
    LOG_TXT += '\nStart merging files:\n'
    # 1. Connect file handles of input and output files for writing
    # Store file handles of input files in a list for iteration
    lst_input_fn = dict_flags['--input']
    lst_input_fh = []  # a list to store file handles of input files
    for fn in lst_input_fn: lst_input_fh.append(xopen(fn, threads=dict_flags['--thread'])) # threads=0 is valid for xopen (ie. no external process is used)
    # File handle for output
    output_fn = dict_flags['--output'] + '.vcf.gz'
    fh_output = PipedCompressionWriter(path=output_fn, threads_flag="-@", program_args=[dict_flags['--write_with']], threads=dict_flags['--thread'])

    # 2. Merge header lines and write to output file (row 0-18 in this version (2021/07,v1))
    # Also Get index number of some columns:
    # - inx_indiv_id_starts: index of the first individual ID (such as R200000348) columns in input .dose.vcf.gz file
    # - inx_info_column: index of INFO columns in input .dose.vcf.gz file
    # - inx_snp_id_column: index of ID columns in input .dose.vcf.gz file
    inx_indiv_id_starts, inx_info_column = merge_header_lines(lst_input_fh, fh_output)
    print('\tMerged header lines (lines start with #)')
    LOG_TXT += '\tMerged header lines (lines start with #)\n'

    # 3. Merge rest lines (from row 19 in this version (2021/07,v1))
    # Merge files, close file handles when done
    merge_files(dict_flags, inx_info_column, inx_indiv_id_starts, lst_input_fh, fh_output)
    for fh in lst_input_fh: fh.close()
    fh_output.close()

    # Print out execution time
    print_execution_time(start_time)

    # Write important processing info into a .log file for user reference
    log_fh = open(dict_flags['--output'] + '.log', 'w')
    log_fh.write(LOG_TXT+'\n')

if __name__ == '__main__': # Called as commandline tool
    run_merge_files('')
# If run as a module in python script, must pass a dictionary arguments to run_merge_files()
