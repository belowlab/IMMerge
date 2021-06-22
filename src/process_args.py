# This function process user input from command line
# Return parsed flags
# Possible flags are:
# --input: (Required) input file names
# --output: (Required) output file names without suffix (eg. BioVU_chr21)
# --verbose: (Optional)
#   If print out detailed information, default is true
#   valid values: [true/True/1, false/False/not 1] v
# --r2_threshold: (Optional)
#   Filtering threshold of imputation quality, default is 0.1.
#   Variants with r2<r2_threshold will be excluded from output
# --r2_output: (Optional)
#   Defines how r2 is calculated in the output. Default is 'first'.
#   Valid values: 'first', 'weighted_average', 'mean'
def process_args(args):
    # For sanity check. Only flags with names in this list are allowed
    lst_flag_names = ['--input', '--output', '--verbose', '--r2_threshold', '--r2_output']

    flag_name = '' # Flag name, such as --input. Initialize with empty string
    flag_val = '' # Flag value, such as chr1.vcf.gz. Initialize with empty string
    lst_flag_vals = [] # A list to store multiple values of --input flag
    if len(args) == 1:  # if no flags provided
        print('Error: No parameters provided, please see --help')
        raise IOError('No parameters provided')
    else:
        dict_flags = dict() # To save flags and flag values
        flag_name = ''
        i = 1 # Skip the first argument
        while i<len(args):
            # The first argument should be a flag (ie. --xxx)
            if i==1 and args[i][0:2]!='--':
                for k, v in dict_flags.items(): print('\t' + k, v)
                raise IOError('Wrong value passed: ' + args[i])
            else: # If not the first one
                if args[i]=='--input': # --input may have multiple values
                    flag_name = args[i]
                    i = i+1
                    while i<len(args) and args[i][0:2]!='--': # Collect everything does not start with '--'
                        lst_flag_vals.append(args[i])
                        i = i+1
                    if len(lst_flag_vals) == 0: # Error if not values for --input
                        for k, v in dict_flags.items(): print('\t' + k, v)
                        raise IOError('Missing values for --input')
                    if len(lst_flag_vals) > 10: # Error if too many values (>10) for --input
                        # May lift this limit
                        for k, v in dict_flags.items(): print('\t' + k, v)
                        raise IOError('To many values for --input. Please merge no more than 10 files at a time')
                    # Assign values to --input flag
                    if dict_flags.get(flag_name) is None:
                        dict_flags[flag_name] = lst_flag_vals
                    else:
                        for k, v in dict_flags.items(): print('\t' + k, v)
                        print('Error: Duplicated flag: '+flag_name)
                        raise IOError('Flag already exist: '+flag_name)

                else: # If not --input flag
                    if args[i][0:2] == '--':
                        flag_name = args[i]
                        if i+1 < len(args): # Make sure there is a value following this flag (--xxx)
                            flag_val = args[i+1]
                            i=i+2
                        else:
                            for k, v in dict_flags.items(): print('\t' + k, v)
                            raise IOError('Missing value for flag: ' + args[i])
                    else: # If not a flag (start with '--'), then this is a wrong input agr
                        for k, v in dict_flags.items(): print('\t' + k, v)
                        raise IOError('Wrong input: '+args[i])

                    # Assign values to --input flag
                    if dict_flags.get(flag_name) is None:
                        dict_flags[flag_name] = flag_val
                    else:
                        for k, v in dict_flags.items(): print('\t' + k, v)
                        print('Error: Duplicated flag: '+flag_name)
                        raise IOError('Flag already exists: '+flag_name)

        # Deal with verbose
        if dict_flags.get('--verbose') is None: dict_flags['--verbose']='true'
        elif dict_flags['--verbose'] in ['true', 'True', '1']: dict_flags['--verbose']='true'
        else: dict_flags['--verbose'] = 'false'

        # Deal with --r2_threshold (default is 0.1)
        if dict_flags.get('--r2_threshold') is None: dict_flags['--r2_threshold']=0.1
        else:
            try:
                dict_flags['--r2_threshold'] = float(dict_flags['--r2_threshold'])
            except:
                print('Error: Invalid value of --r2_threshold:', dict_flags['--r2_threshold'])
                print('\tValue of --r2_threshold should be numeric and between 0 and 1')
                raise IOError('Invalid value of --r2_threshold')

            if dict_flags['--r2_threshold'] < 0 or dict_flags['--r2_threshold'] > 1:
                print('Error: Invalid value of --r2_threshold:', dict_flags['--r2_threshold'])
                print('\tValue of --r2_threshold should be numeric between 0 and 1')
                raise IOError('Invalid value of --r2_threshold')

        # Deal with --r2_output
        if dict_flags.get('--r2_output') is None: dict_flags['--r2_output']='first'
        else:
            if dict_flags['--r2_output'] not in ['first', 'weighted_average', 'mean']:
                print('Error: Invalid value of --r2_output:', dict_flags['--r2_output'])
                print('Value of --r2_output should be: first, weighted_average or mean')
                raise IOError('Invalid value of --r2_output')

        # If no error raised up to here, print flags and values used
        print('\nInput options used:')
        for k, v in dict_flags.items():
            if k not in lst_flag_names:
                print('Error: Unrecognized flag: '+k)
                print('Only these flags are allowed:', lst_flag_names)
                raise IOError('Invalid flag: ' + flag_name)
            else: print('\t'+k, v)
        return dict_flags

