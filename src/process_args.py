# This function process user input from command line
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
def process_args(args):
    # For sanity check. Only flags with names in this list are allowed
    lst_flag_names = ['--input', '--output', '--thread', '--missing', '--na_rep',
                      '--r2_threshold', '--r2_output', '--duplicate_id', '--help']
    # Save flags and flag values in a dictionary
    dict_flags = dict()

    # flag_name = '' # Flag name, such as --input. Initialize with empty string
    # flag_val = '' # Flag value, such as chr1.vcf.gz. Initialize with empty string
    lst_flag_vals = [] # A list to store multiple values of --input flag
    if len(args) == 1:  # if no flags provided
        print('Error: No parameters provided, please see --help\nExit')
        exit()
    else:
        i = 1 # Skip the first argument
        while i<len(args):
            # The first argument should be a flag (ie. --xxx)
            if i==1 and args[i][0:2]!='--':
                for k, v in dict_flags.items(): print('\t' + k, v)
                raise IOError('Wrong value passed: ' + args[i])
            else: # If not the first one
                if args[i] == '--help': # Print out help info, ignore other parameters (if provided)
                    import print_help
                elif args[i]=='--input': # --input may have multiple values
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
                        print('Error: Duplicated flag: '+flag_name, '\n')
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
                        print('Error: Duplicated flag: '+flag_name, '\n')
                        raise IOError('Flag already exists: '+flag_name)

        # Check --input
        if dict_flags.get('--input') is None:
            print('Error: Missing input files\n')
            raise IOError('Required flag missing: --input')

        # Check --output
        if dict_flags.get('--output') is None: dict_flags['--output']='merged'

        # Check --thread
        if dict_flags.get('--thread') is None: dict_flags['--thread']=1
        else:
            try:
                dict_flags['--thread'] = int(dict_flags['--thread'])
            except:
                print('Error: Invalid value of --thread:', dict_flags['--thread'])
                print('\tValue of --thread should be an integer\n')
                raise IOError('Invalid value of --thread')
            # Assign 1 to --thread if user supplied a value<0
            if dict_flags['--thread']<0: dict_flags['--thread']=1

        # Check --missing
        if dict_flags.get('--missing') is None: dict_flags['--missing'] = 0
        else:
            try:
                dict_flags['--missing'] = int(dict_flags['--missing'])
            except:
                print('Error: Invalid value of --missing:', dict_flags['--missing'])
                print('\tValue of --missing should be an integer\n')
                raise IOError('Invalid value of -cat ../ou  -missing')
            if dict_flags['--missing'] < 0 or dict_flags['--missing'] > len(dict_flags['--input']):
                print('Error: Invalid value of --missing:', dict_flags['--missing'])
                print('\tValue of --missing should be an integer between 0 and number of input files\n')
                raise IOError('Invalid value of --missing')

        # Check --na_rep
        # Default is . as VCFs uses "." to represent missing values, trailing fields after missing genotyoe can be ignored
        # In the output VCF, missing genotype field is ".|.", other fileds are "."
        if dict_flags.get('--na_rep') is None: dict_flags['--na_rep'] = '.'

        # Check --r2_threshold (default is 0)
        if dict_flags.get('--r2_threshold') is None: dict_flags['--r2_threshold'] = 0
        else:
            try:
                dict_flags['--r2_threshold'] = float(dict_flags['--r2_threshold'])
            except:
                print('Error: Invalid value of --r2_threshold:', dict_flags['--r2_threshold'])
                print('\tValue of --r2_threshold should be numeric and between 0 and 1\n')
                raise IOError('Invalid value of --r2_threshold')

            if dict_flags['--r2_threshold'] < 0 or dict_flags['--r2_threshold'] > 1:
                print('Error: Invalid value of --r2_threshold:', dict_flags['--r2_threshold'])
                print('\tValue of --r2_threshold should be numeric between 0 and 1\n')
                raise IOError('Invalid value of --r2_threshold')

        # Check --r2_output
        if dict_flags.get('--r2_output') is None: dict_flags['--r2_output']='first'
        else:
            if dict_flags['--r2_output'] not in ['first', 'weighted_average', 'mean', 'min', 'max']:
                print('Error: Invalid value of --r2_output:', dict_flags['--r2_output'])
                print('Value of --r2_output should be: first, weighted_average or mean\n')
                raise IOError('Invalid value of --r2_output')

        # Check --duplicate_id
        if dict_flags.get('--duplicate_id') is None: dict_flags['--duplicate_id'] = 0
        else:
            try:
                dict_flags['--duplicate_id'] = int(dict_flags['--duplicate_id'])
            except:
                print('Error: Invalid value of --duplicate_id:', dict_flags['--duplicate_id'])
                print('\tValue of --duplicate_id should be ann integer\n')
                raise IOError('Invalid value of --duplicate_id')
            if dict_flags['--duplicate_id'] < 0:
                print('Error: Invalid value of --duplicate_id:', dict_flags['--duplicate_id'])
                print('\tValue of --duplicate_id should be an integer between 0 and number of individuals\n')
                raise IOError('Invalid value of --duplicate_id')

        # If no error raised up to here, print flags and values used
        print('\nInput options used:')
        for k, v in dict_flags.items():
            if k not in lst_flag_names:
                print('Error: Unrecognized flag: '+k)
                print('Only these flags are allowed:', lst_flag_names, '\n')
                raise IOError('Invalid flag: ' + k)
            else:
                if k=='--na_rep' and dict_flags['--missing']==0:
                    print('\t' + k, v, '(ignored since --missing is 0)')
                elif k=='--input' and len(dict_flags['--input'])<2:
                    print('Error: At least two files are needed for --input:', dict_flags['--input'])
                    raise IOError('Missing --input files')
                else: print('\t'+k, v)
        return dict_flags

