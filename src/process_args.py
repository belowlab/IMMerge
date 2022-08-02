# This script process user input from command line
# Return parsed flags in a dictionary (dict_flags)
# Possible flags are:
# --input: (Required) input file names
# --output: (Required) output file names without suffix (eg. BioVU_chr21)
# --thread: (Optional)
#   Default is 1
#   Defines how many thread to use in multiprocessing
#   valid values: int >= 1
# --missing: (optional, int)
#   Defines number of missing values allowed for each variant.
#   Default is 0. Cannot exceed total number of input files.
#   If --missing is 0, only variants shared by all input files will be saved in merged result
# --na_rep: (optional)
#   Defines what symbol to use for missing values. Default is 'NA'
#   This flag will be ignored if --missing is 0.
# --r2_threshold: (Optional)
#   Filtering threshold of imputation quality, default is 0 (no filtering).
#   Variants with r2<r2_threshold will be excluded from output
# --r2_output: (Optional)
#   Defines how r2 is calculated in the output. Default is 'first'.
#   Valid values: 'first', 'weighted_average', 'mean'
#   If use 'first', --missing must be set to 0
# --duplicate_id: (optional, int)
#   Default is 0. Defines number of duplicated individuals in each input file.
#   Duplicated IDs should be the first N columns in each file and not mixed with unique IDs.
#   Starting from the second input file, data of the first N individuals will be skipped in the merged output
# --write_with: (optional)
#   Default is bcftools. Write to bgziped file with bcftools. User can supply specific path to bcftools such as /user/bin/bcftools
#   Use --write_with "/data100t1/gapps/htslib-1.9/bgzip" for our server
import argparse
import os

def process_args():
    parser = argparse.ArgumentParser()
    lst_args = ['--input', '--info', '--output', '--thread', '--missing', '--na_rep', '--r2_threshold', '--r2_output',
                '--r2_cap', '--duplicate_id', '--check_duplicate_id', '--write_with', '--meta_info', '--verbose']
    # Help messages of each option
    dict_help = {
        '--input': '(Required) Files to be merged, multiple files are allowed. Must in gzipped or bgziped VCF format',
        '--info': '(Optional) Directories and name to info files. Default path is the same directory as corresponding input file, default info file share the same name as input file, except for suffix (.info.gz)',
        '--output': '(Optional) Default is "merged". output file name without suffix',
        '--thread': '(Optional) Default value is 1. Defines how many thread to use in multiprocessing. If number of threads <0, will use 1 instead of user supplied value.',
        '--missing': '(Optional) Default is 0. Defines number of missing values allowed for each variant',
        '--na_rep': '(Optional) Default is "." (ie. ".|." for genotype values). Defines what symbol to use for missing values. This flag is ignored if --missing is 0',
        '--r2_threshold': '(Optional) Default is 0, ie. no filtering. Only variants with combined imputation quality r2≥r2_threshold will be saved in the merged file',
        '--r2_output': '(Optional) Default is "z_transformation". Defines how imputation quality score is calculated in the output file.',
        '--r2_cap': '(Optional) Adjust R squared by --r2_cap if imputation quality Rsq=1. Only valid for z transformation to avoid infinity',
        '--duplicate_id': '(Optional) Default is 0. Defines number of duplicated individuals in each input file. Duplicated IDs should be the first N columns in each file',
        '--check_duplicate_id':'(Optional) Default is False. Check if there are duplicated IDs, then rename non-first IDs to ID:2, ID:3, ..., ID:index_of_input_file+1',
        '--write_with': '(Optional) Default is bgzip. Write to bgziped file with bgzip. User can supply specific path to bgzip such as /user/bin/bgzip',
        '--meta_info': "(Optional) Valid values are {index of input file (1-based), 'none', 'all'}. What meta information (lines start with '##') to include in output file. Default is 1 (meta information from the first input file)",
        '--verbose': '(Optional) Default is False. Print more messages.'}

    # Default values and data types of optional flags
    # --r2_output and --verbose also have choices from limited values
    dict_default = {'--info': ['', str],
                    '--output': ['merged', str],
                    '--thread': [1, int],
                    '--missing': [0, int],
                    '--na_rep': ['.', str],
                    '--r2_threshold': [0, float],
                    '--r2_output': ['z_transformation', str, ['first', 'mean', 'weighted_average', 'z_transformation', 'min', 'max']],
                    '--r2_cap': [10e-4, float],
                    '--check_duplicate_id':['0', str, ['false', '0', 'False', 'true', 'True', '1']],
                    '--duplicate_id': [0, int],
                    '--write_with':['bgzip', str],
                    '--meta_info':['1', str],
                    '--verbose': ['0', str, ['false', '0', 'False', 'true', 'True', '1']]}
    # Add arguments
    for arg in lst_args:  # If user provide arguments not in the list, they will not be used (and no error message)
        if arg == '--input':
            parser.add_argument(arg, help=dict_help[arg], nargs='*', required=True)
        elif arg == '--info':
            parser.add_argument(arg, help=dict_help[arg], nargs='*', default='')
        elif arg == '--r2_output' or arg == '--verbose' or arg=='--check_duplicate_id':
            parser.add_argument(arg, help=dict_help[arg], default=dict_default[arg][0], type=dict_default[arg][1], choices=dict_default[arg][2])
        else:
            parser.add_argument(arg, help=dict_help[arg], default=dict_default[arg][0], type=dict_default[arg][1])

    args = parser.parse_args()

    # Convert --check_duplicate_id to boolean
    if args.check_duplicate_id.upper()=='FALSE' or args.check_duplicate_id=='0':
        args.check_duplicate_id = False
    else:
        args.check_duplicate_id = True

    # Convert verbose to boolean
    if args.verbose.upper()=='FALSE' or args.verbose=='0':
        args.verbose = False
    else:
        args.verbose = True

    # Add values to --info if not provided by user
    if args.info=='':
        args.info = []
        for val in args.input:
            path, fn = os.path.split(val)
            info_fn = fn.split('.')[0] + '.info.gz'
            args.info.append(os.path.join(path, info_fn))

    dict_flags = {}  # Store arguments in a dictionary to return
    for arg in lst_args:
        dict_flags[arg] = eval('args.' + arg[2:])

    # Print flags and values loaded
    for k, v in dict_flags.items():
        if k == '--na_rep' and dict_flags['--missing'] == 0:
            print('\t' + k, v, '(ignored since --missing is 0)')
        elif k == '--r2_cap' and dict_flags['--r2_output'] != 'z_transformation':
            print('\t' + k, v, '(ignored since --r2_output is not z_transformation)')
        else:
            print('\t' + k, v)

    # ################ Sanity checks (type check is taken care of by type argument of add_argument) ################
    # Check --input
    if len(dict_flags['--input']) <= 1:  # Number of input files
        print('Error: Invalid value of --input:', dict_flags['--input'])
        print('\t- At least two files are needed to merge\nExit')
        exit()
    for val in dict_flags['--input']: # Check if file exists
        if not os.path.isfile(val):
            print(f"Error: Input file not found: {val}")
            exit()

    # Check --info
    if len(dict_flags['--input']) != len(dict_flags['--info']):
        print(f"Error: number of input files ({len(dict_flags['--input'])}) does not match number of info files ({len(dict_flags['--info'])}).\nExit")
        exit()
    for info_file in dict_flags['--info']:
        if not os.path.isfile(info_file):
            print(f"Error: Info file ({info_file}) does not exist.\nExit")
            exit()

    # Check --output
    path, fn = os.path.split(dict_flags['--output'])
    if not os.path.isdir(path): # Create output directory if does not exist
        os.mkdir(path)

    # Check --thread
    if dict_flags['--thread'] <= 0: dict_flags['--thread'] = 1  # Assign 1 to --thread if user supplied a value<=0

    # Check --missing
    if dict_flags['--missing'] < 0 or dict_flags['--missing'] > len(dict_flags['--input']):
        print('Error: Invalid value of --missing:', dict_flags['--missing'])
        print('\t- Value of --missing should be an integer between 0 and number of input files\nExit')
        exit()

    # Check --r2_threshold (default is 0)
    if dict_flags['--r2_threshold'] < 0 or dict_flags['--r2_threshold'] > 1:
        print('Error: Invalid value of --r2_threshold:', dict_flags['--r2_threshold'])
        print('\t- Value of --r2_threshold should be numeric between 0 and 1\nExit')
        exit()

    # # Check --r2_output
    # if dict_flags['--r2_output'] not in ['first', 'weighted_average', 'z_transformation', 'mean', 'min', 'max']:
    #     print('Error: Invalid value of --r2_output:', dict_flags['--r2_output'])
    #     print('\t- Value of --r2_output should be: first, weighted_average, z_transformation, mean, min or max\nExit')
    #     exit()

    # Check --r2_cap
    if dict_flags['--r2_cap'] <0 or dict_flags['--r2_cap']>=1:
        print('Error: Invalid value of --r2_cap:', dict_flags['--r2_cap'])
        print('\t- Value of --r2_cap should be numeric between 0 and 1\nExit')
        exit()

    # Check --duplicate_id
    if dict_flags['--duplicate_id'] < 0:
        print('Error: Invalid value of --duplicate_id:', dict_flags['--duplicate_id'])
        print('\t- Value of --duplicate_id should be an integer between 0 and number of individuals\nExit')
        exit()

    # Check --meta_info
    meta_flag = True
    if dict_flags['--meta_info'].isnumeric():
        dict_flags['--meta_info'] = int(dict_flags['--meta_info'])
        if dict_flags['--meta_info'] <= 0 or dict_flags['--meta_info']>len(dict_flags['--input'] ):
            meta_flag = False
    elif dict_flags['--meta_info'] != 'none' and dict_flags['--meta_info'] != 'all':
        meta_flag = False
    if not meta_flag:
        print('Error: Invalid value of --meta_info:', dict_flags['--meta_info'])
        print('\t- Value of --meta_info should be an integer indicating index (1-based) of input file, "none", or "all"\nExit')
        exit()

    # If no error raised up to here, return flags
    return dict_flags

