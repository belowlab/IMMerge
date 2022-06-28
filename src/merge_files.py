# Step 02
# This code is used to merge imputed files from TopMed based on list of saved variants in prior steps
# The goal is to combine post-imputation vcf.gz files into one vcf.gz file

from xopen import PipedCompressionWriter, xopen # Faster than gzip
# import sys
import process_args
import get_SNP_list
import check_r2_setting_for_imputation
import time
import multiprocessing
import bgzip # vcf.gz file is actually bgzip file not gzip file, but can be opened with gzip

print('\nVersion 1.0')
start_time = time.time()    # Track execution time
print('Job starts at (local time) ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))

# File name can be passed to this code in terminal, or use import this code as in a script (need to modify a little)
# args = sys.argv

# Process terminal input
# dict_flags contains values for below flags:
#   --input, --output, --thread, --missing, --r2_threshold, --r2_output
# dict_flags = process_args_ongoing.process_args(args) # Process terminal input
dict_flags = process_args.process_args() # Process terminal input

# Write some info into a log file
log_fn = dict_flags['--output'] + '.log' # Save Important processing info into a .log file for user reference
with open(log_fn, 'w') as log_fh:
    # Write start time
    log_fh.write('Job starts at (local time) ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)) + '\n\n')
    log_fh.write('Input options used:\n')
    for k,v in dict_flags.items():
        log_fh.write('\t'+ k+': ' + str(v) + '\n')

    # write info.gz file names
    log_fh.write('\nBelow .info.gz files were used:\n')
    for fn in dict_flags['--input']:
        info_fn = fn.split('/')[-1].split('.')[0] + '.info.gz'
        log_fh.write('\t' + info_fn + '\n')

    # Write number of variants and individuals to be processed
    log_fh.write('\nNumber of variants:\n') # Actual numbers are written in get_SNP_list.py
    # log_fh.write('\tTotal number of all input files combined:')
    # log_fh.write('\tNumber of saved variants:')
    # log_fh.write('\tNumber of excluded variants:')

check_r2_setting_for_imputation.check_imputatilson_parameters(lst_fn=dict_flags['--input'])
lst_number_of_individuals = get_SNP_list.get_snp_list(dict_flags)

# ----------------------- Helper functions -----------------------
# This function prints execution time at the end of run
def print_execution_time(satrt_time):
    # merge_files()
    duration = time.time() - start_time
    hours = duration // 3600
    duration = duration % 3600
    minutes = duration // 60
    seconds = duration % 60
    print('\nDuration: {:.2f} hours, {:.2f} minutes, {:.2f} seconds'.format(hours, minutes, seconds))
    print('End of run')

    # Write into log file
    log_fn = dict_flags['--output'] + '.log'
    with open(log_fn, 'a') as log_fh:
        log_fh.write('\n\nDuration: {:.2f} hours, {:.2f} minutes, {:.2f} seconds'.format(hours, minutes, seconds))
        log_fh.write('\n\nEnd of run')

# Merge header lines of input files (.dose.vcf.gz) together and write into a BGZF file (bgzip file)
# Parameters:
# - number_of_dup_id：User provided number of duplicated IDs in each sub input file
# - fh_output_raw: file handle of output file (need to use bgzip module for actual writing)
# - lst_input_fh: list of file handles of input files
def merge_header_lines(lst_input_fh, fh_output, number_of_dup_id=dict_flags['--duplicate_id']):
    # # Read in SNPs (with corresponding combined MAF and r2) need to be kept from file
    # lst_snp, lst_alt_frq, lst_maf, lst_r2 = get_snp_and_info_lst()

    # Read trough header lines (line start with ##) and write into output file
    # Process all files at the same time. All input files should have the same number of header lines
    line = ''
    for fh in lst_input_fh:  # Read the 1st line of all input files
        line = fh.readline().strip()

    # Use encode() to convert string to byte
    # This is to make sure bgzip module works, since it uses bytearray to write
    fh_output.write(line)
    line = lst_input_fh[0].readline().strip()  # Read the 2nd line of the first input files

    inx_indiv_id_starts = 0  # From which column genotype data starts
    inx_info_column = 0 # Find index of column "INFO", in order to replace values later
    while line[0] == '#': # Info and header lines all start with "#"
        if line[0:2] == '##':  # If Info lines, write everything directly
            fh_output.write(('\n' + line))
            for fh in lst_input_fh[1:]:  # Read the rest file handles, but do not need to write them
                line = fh.readline().strip()
        else:  # If column header line, then merge and write (line starts with single #)
            # Add a few lines about how this file is generated (with '##')
            extra_info_line = '\n##VCFmerge=<Input files='
            for i in dict_flags['--input']: # Write input file names
                if '/' in i: # remove directory if possible, only use file names
                    i = i.split('/')[-1]
                extra_info_line = extra_info_line + i + ', '
            extra_info_line = extra_info_line[:-2]+'>' # Remove the last ', '
            fh_output.write((extra_info_line + '\n'))
            fh_output.write(('##VCFmerge=<Date=' + time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime()) + '>\n'))
            fh_output.write(('##VCFmerge=<Missing=' + str(dict_flags['--missing']) + ', duplicate_id='+str(dict_flags['--duplicate_id'])+', na_rep='+dict_flags['--na_rep']+'>\n'))
            # fh_output.write(('##VCFmerge=<Rsq filter=' + str(dict_flags['--r2_threshold']) + '>\n'))
            fh_output.write(('##VCFmerge=<Merged MAF=weighted average, Merged AF=weighted average>\n'))
            fh_output.write(('##VCFmerge=<Merged Rsq=' + dict_flags['--r2_output'] + ', Rsq filter=' + str(
                dict_flags['--r2_threshold']) + '>\n+'))

            # In current TOPMed version these are columns of shared information, then followed by individual IDs:
            #   - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
            inx_indiv_id_starts = line.split().index('FORMAT') + 1  # Find from which column genotype data starts
            inx_info_column = line.split().index('INFO')
            # inx_snp_id_column = line.split().index('ID')
            print('\tIndividual genotype data starts from column #:', inx_indiv_id_starts + 1)
            # Create merged header line, initialize with header from the first file
            line_to_write = line
            for fh in lst_input_fh[1:]: # Read header line from rest file handles
                line = fh.readline().strip()
                # Set maxsplit in slit() function to remove duplicated IDs
                line_to_write = line_to_write + '\t' + line.split(maxsplit=number_of_dup_id + inx_indiv_id_starts)[-1]
            fh_output.write((line_to_write + '\n'))
            break

        line = lst_input_fh[0].readline().strip()  # Read line of the first file to decide what to do

    return inx_indiv_id_starts, inx_info_column

# This function read line of a given file handle, until find line of given SNP (by parameter snp))
# Return string of that line (if not found then return an empty string '')
def search_SNP_and_read_lines(snp, fh):
    inx_snp_ID = 2  # Index of column that contains variant ID such as chr21:10000777:C:A (value is 2 in 2021/07 version)
    line = fh.readline().strip()  # Read in each line to search for the given variant
    # input_snp = line.split(maxsplit=inx_snp_ID + 1)[-2]
    while line != '':  # Keep read line until find given snp
        input_snp = line.split(maxsplit=inx_snp_ID + 1)[-2]
        if input_snp == snp: # if found
            break
        else:
            line = fh.readline().strip()  # Read in each line to search for the given variant
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
    # Load variants_kept.txt as a reference of merged file
    # Then compare snp ids in each input file with the SNPs in variants_kept.txt
    # And replace ALT_frq, MAF and r2 values in INFO column with new values calculated from all input files
    # - ie, use these columns in variant_kept.txt file: ALT_Frq_combined, MAF_combined, Rsq_combined, Genotyped
    count = 0  # For console output
    fh_snp_kept = open(dict_flags['--output']+'_variants_kept.txt')
    print('\t'+dict_flags['--output']+'_variants_kept.txt file loaded')

    # Read the first line (header), find index of combined ALT_freq, MAF, Rsq
    line = fh_snp_kept.readline().strip().split()
    indx_ALT_Frq_combined = line.index('ALT_Frq_combined')
    indx_MAF_combined = line.index('MAF_combined')
    indx_Rsq_combined = line.index('Rsq_combined')

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
            genotyped = lst_snp_kept[3] # "Genotyped" is always the 4th column, so can be hard coded
            # Create new INFO value to replace INFO column in merged line
            new_info_val = 'AF='+ALT_Frq_combined+';MAF='+MAF_combined+';R2='+Rsq_combined+';'+genotyped

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
            if count%5000 == 0: print('. {} variants merged'.format(count))
            elif count%100 == 0: print('.', end='', flush=True)
            elif count%5000 == 1: print('\t', end='', flush=True)
        else:
            print('\n\tEnd of snp to be kept file, total of {} variants merged'.format(count))
            break  # If reach end of SNP to be kept, then stop
    # Close file handles when done
    fh_snp_kept.close()

# ----------------------- End of helper functions -----------------------

# A function to run each step all together
def run_merge_files():
    print('\nStart merging files:')
    # 1. Connect file handles of input and output files for writing
    # Store file handles of input files in a list for iteration
    lst_input_fn = dict_flags['--input']
    lst_input_fh = []  # a list to store file handles of input files
    for fn in lst_input_fn: lst_input_fh.append(xopen(fn, threads=dict_flags['--thread']-1)) # threads=0 is valid for xopen (ie. no external process is used)
    # File handle for output
    # output_fn = dict_flags['--output'] + '.dose.vcf.gz' # Do not always need this dose suffix
    output_fn = dict_flags['--output'] + '.vcf.gz'
    fh_output = PipedCompressionWriter(path=output_fn, threads_flag="-@", program_args=["bgzip", f"-I {output_fn}i"], threads=dict_flags['--thread'])

    # 2. Merge header lines and write to output file (row 0-18 in this version (2021/07,v1))
    # Also Get index number of some columns:
    # - inx_indiv_id_starts: index of the first individual ID (such as R200000348) columns in input .dose.vcf.gz file
    # - inx_info_column: index of INFO columns in input .dose.vcf.gz file
    # - inx_snp_id_column: index of ID columns in input .dose.vcf.gz file
    inx_indiv_id_starts, inx_info_column = merge_header_lines(lst_input_fh, fh_output)
    print('\tMerged header lines (lines start with #)')

    # 3. Merge rest lines (from row 19 in this version (2021/07,v1))
    # Merge files, close file handles when done
    merge_files(dict_flags, inx_info_column, inx_indiv_id_starts, lst_input_fh, fh_output)
    for fh in lst_input_fh: fh.close()
    fh_output.close()

    # Print out execution time
    print_execution_time(start_time)

# ----------------- Profile memory usage ---------------------
'''from memory_profiler import profile

@profile
def run_merge_script():
    multiprocessing.set_start_method("fork")  # This is necessary for python 3.8, but won't matter in other versions
    # Use user defined number of cores to do multi processing, unless not defined
    if dict_flags['--thread'] == 1:  # no multiprocessing
        run_merge_files()
    elif dict_flags['--thread'] <= multiprocessing.cpu_count():
        number_of_cores_to_use = dict_flags['--thread']
        with multiprocessing.Pool(number_of_cores_to_use) as p:
            run_merge_files()
    else:
        number_of_cores_to_use = multiprocessing.cpu_count()
        with multiprocessing.Pool(number_of_cores_to_use) as p:
            run_merge_files()

if __name__ == '__main__':
    run_merge_script()'''

if __name__ == '__main__':
    run_merge_files()
