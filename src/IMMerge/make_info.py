# This script create .info.gz file from bgzipped vcf file.
# info files are already included if all output files are downloaded from TOPmed imputation server
# Requires bcftools
import argparse
import os
from xopen import PipedCompressionWriter, xopen

def process_args(arg_list = ''):
    '''
    Process arguments from command line or passed as a list from python module call.
    Return args object that contains all arguments and corresponding values
    '''

    lst_arg = ['--input', '--output_dir', '--output_fn', '--col_names', '--thread', '--write_with', '--use_rsid',
               '--mixed_genotype_status', '--genotyped', '--imputed', '--verbose'] # Available arguments
    parser = argparse.ArgumentParser(description='Generate info files')
    parser.add_argument('--input', nargs='*', type=str, required=True,
                        help='(Required) Multiple input files are allowed. Must in gzipped or bgziped VCF format.')
    parser.add_argument('--output_dir', type=str, default='',
                        help='(Optional) Directory for output files. Default is current working directory.')
    parser.add_argument('--output_fn', type=str, nargs='*',
                        help='(Optional) Default is input file name with suffix replaced by ".info.gz".')
    parser.add_argument('--col_names', nargs='*', default=['AF', 'MAF', 'R2', 'IMPUTED/TYPED/TYPED_ONLY'],
                        help="(Optional) Default is ['AF', 'MAF', 'R2', 'IMPUTED/TYPED/TYPED_ONLY']. Column names of alt frequency, MAF, imputation quality score, genotype status. Separated by space.")
    parser.add_argument('--thread', default=1, type=int,
                        help='(Optional) Default is 1. Defines how many thread to use in multiprocessing. If number of threads <0, will use 1 instead of user supplied value.')
    parser.add_argument('--write_with', default='bgzip', type=str,
                        help='(Optional) User specified program to write compressed file. Default is bgzip. (gzip can be used unless required to bgzip info files)')
    parser.add_argument('--use_rsid', default='False', choices=['false', '0', 'False', 'true', 'True', '1'], type=str,
                        help='(Optional) Default is False. If input VCFs use rsID instead of chr:pos:ref:alt, set this option to True to avoid duplicate IDs (rsID may not be unique). Also need to use the same setting in merging step.')
    parser.add_argument('--mixed_genotype_status', default='False', type=str, choices=['false', '0', 'False', 'true', 'True', '1'],
                        help='(Optional) Default is False. Valid values are (not case-sensitive): {0|1|True|False}. Whether some variants have more than one genotype status (True) or not (False). Use together with arguments --genotyped and --imputed. If False then output genotype status of each variant is the last genotype status in its INFO column. If True then output genotype status fo reach variant will be: ALL=all genotyped, SOME=at least one genotyped, NONE=no genotyped.')
    parser.add_argument('--genotyped', default='TYPED/TYPED_ONLY', type=str,
                        help='(Optional) Default is TYPED/TYPED_ONLY in concordance to TOPMed output. Label for genotyped variants. Multiple values can be supplied in one string separated by /. Only evaluated when --mixed_genotype_status is True.')
    parser.add_argument('--imputed', default='IMPUTED', type=str,
                        help='(Optional) Default is IMPUTED in concordance to TOPMed output. Label for imputed variants. Multiple values can be supplied in one string separated by /. Only evaluated when --mixed_genotype_status is True.')
    parser.add_argument('--verbose', default='False', choices=['false', '0', 'False', 'true', 'True', '1'], type=str,
                        help='(Optional) Print help message if True. Default is False. Valid values are (not case-sensitive): {0|1|True|False}.')

    if arg_list == '':
        args = parser.parse_args() # if it is called from command line
    else: args = parser.parse_args(arg_list) # If it is used as a python module

    # Turn --verbose to boolean
    if args.verbose.upper()=='FALSE' or args.verbose=='0':
        args.verbose = False
    else:
        args.verbose = True
        parser.print_help()  # Print help message if desired

    # Turn --use_rsid to boolean
    if args.use_rsid.upper()=='FALSE' or args.use_rsid=='0':
        args.use_rsid = False
    else:
        args.use_rsid = True

    # Turn --mixed_genotype_status to boolean
    if args.mixed_genotype_status.upper()=='FALSE' or args.mixed_genotype_status=='0':
        args.mixed_genotype_status = False
    else:
        args.mixed_genotype_status = True

    if args.output_fn is None: # Fill output file names as input name + '.info.gz'
        args.output_fn = [os.path.split(fn)[1].split('.')[0] + '.info.gz' for fn in args.input]
    if args.thread<0: args.thread = 1

    print('\nOptions used')
    for arg in lst_arg:
        print('\t' + arg, eval('args.'+arg[2:]))
    print('\n')

    # ################################ Sanity checks ################################
    if args.verbose:
        print('#### Sanity check: checking user inputs')

    invalid_flag = False # Track if any input file is invalid
    # Check if input files exist
    for fn in args.input:
        if not os.path.isfile(fn):
            print(f'Error: Invalid input file name: {fn} is not a file or does not exist')
            invalid_flag = True

    # Check if number of input and output files matches
    if len(args.input) != len(args.output_fn):
        print(f'Error: Number of input file ({len(args.input)}) does not match number of output file ({len(args.output_fn)})')
        invalid_flag = True

    if len(args.col_names) > 4:
        print(f'Error: Too many values. Only need (4) column names of alt frequency, MAF, imputation quality score and genotyped.')
        invalid_flag = True
    if len(args.col_names) < 4:
        print(f'Error: Not enough values. Need (4) column names of alt frequency, MAF, imputation quality score and genotyped. Values should be separated by space.')
        invalid_flag = True

    if invalid_flag:
        print('Exit')
        exit()

    if args.verbose:
        print('#### Sanity check: Done. All inputs are valid')
    # ################################ End of sanity check ################################
    return args

def make_args_list(args_dict):
    '''
    Make a list of arguments to pass to process_args when IMMerge is used as a python module
    Params:
    - args_dict: a dictionary of arguments similar to {'--input':[input1, input2]}
    Return:
    A list of arguments (to pass to parse_args())
    '''
    args_list = []
    for k, v in args_dict.items():
        args_list.append(k)
        if isinstance(v, list):
            for val in v: args_list.append(val)
        else: args_list.append(v)

    return args_list

def genotype_status(gp_list, genotyped_label, imputed_label):
    '''
    Handle mixed genotype status. Return final genotype status based a list of genotype status
    Params:
    - gp_list: a list of genotype status such as ['IMPUTED', 'IMPUTED', 'GENOTYPED']
    - genotyped_label: label(s) for genotype status "genotyped", separated by '/'
    - imputed_label: label(s) for genotype status "imputed", separated by '/'
    Return:
    A string or -1. ALL=all genotyped, SOME=at least one genotyped, NONE=no genotyped, -1=field missing
    '''
    genotyped_list = genotyped_label.split('/')
    imputed_list = imputed_label.split('/')

    # Crate flags to track if there are certain types of genotype status
    status, any_gped, any_imputed = -1, False, False
    for gp_status in gp_list:
        if gp_status in genotyped_list:
            any_gped = True
        elif gp_status in imputed_list:
            any_imputed = True
    if any_gped and any_imputed:
        status = 'SOME'
    elif any_gped and not any_imputed:
        status = 'ALL'
    elif not any_gped and any_imputed:
        status = 'NONE'
    return status

def write_info(args):
    '''
    Write information to output file (.info.gz)
    Params:
    - args: a dictionary of arguments as {'--input':['file1.vcf.gz', 'file2.vcf.gz'], '--output_fn':'output.info.gz'}. Valid flags are:
        --input: (Required) Multiple input files are allowed. Must in gzipped or bgziped VCF format.
        --output_dir: (Optional) Directory for output files. Default is current working directory.
        --output_fn: (Optional) Default is input file name with suffix replaced by ".info.gz".
        --col_names: (Optional) Default is ['AF', 'MAF', 'R2', 'IMPUTED/TYPED/TYPED_ONLY'].
                    Column names of alt frequency, MAF, imputation quality score, genotyped. Separated by space.
        --thread: (Optional) Default is 1. Defines how many thread to use in multiprocessing.
                If number of threads <0, will use 1 instead of user supplied value.')
        --write_with: (Optional) Default is bgzip, but gzip is also valid. User specified program to write compressed file.
        --use_rsid: (Optional) Default is False.
                    If input VCFs use rsID instead of chr:pos:ref:alt, set this option to True to avoid duplicate IDs (rsID may not be unique).
                    Also need to use the same setting in merging step.
        --verbose: (Optional) Print help message if True. Default is False. Valid values are (not case-sensitive): {0|1|True|False}

    Return:
    - A .info.gz file is created as settings used in '--output_dir' and '--output_fn'
    '''
    if args == '': # If called from command line
        args = process_args()
    else: # If called as a module
        args_list = make_args_list(args) # Convert the dictionary of arguments to a list to be parsed
        args = process_args(args_list)

    # Create directory if user provide --output_dir does not exist
    if len(args.output_dir) and (not os.path.isdir(args.output_dir)):
        if args.verbose:
            print(f'#### Output directory {args.output_dir} does not exist. Creating output directory')
        os.mkdir(args.output_dir)

    for i in range(len(args.input)):
        if args.verbose:
            print(
                f'\n#### Processing files: Start. input file = {args.input[i]}; output file = {os.path.join(args.output_dir, args.output_fn[i])}')

        output_fh = PipedCompressionWriter(path=os.path.join(args.output_dir, args.output_fn[i]), threads_flag="-@",
                                           program_args=[args.write_with],
                                           threads=args.thread)
        # Need these columns: 'SNP', 'REF(0)', 'ALT(1)', 'Genotyped','ALT_Frq', 'MAF', 'Rsq'
        output_fh.write('SNP\trsID\tREF(0)\tALT(1)\tALT_Frq\tMAF\tRsq\tGenotyped\n')
        if args.verbose:
            print(f'#### Processing files: Write column headers to output file {args.output_fn[i]}')

        with xopen(args.input[i]) as input_fh:
            line = input_fh.readline().strip()
            while line[:2] == '##':  # Read through header lines start with ##
                line = input_fh.readline().strip()
            column_names = line.split(maxsplit=10)  # The header line has 8 fixed columns, so do not need to split the entire line
            if args.use_rsid:  # If ID in input VCF is rsID, need to extract chr, pos, ref and alt instead of ID column
                # Find index of #CHROM, POS, REF and ALT
                idx_chr = column_names.index('#CHROM')
                idx_pos = column_names.index('POS')
                idx_rsid = column_names.index('ID')  # In this case ID column contains rsID

            # Find index of other desired columns
            idx_snp = column_names.index('ID')
            idx_ref = column_names.index('REF')
            idx_alt = column_names.index('ALT')
            idx_info = column_names.index('INFO')  # INFO column

            # Get indices of AF (alt_frq), MAF, Rsq and genotyped from INFO column
            line = input_fh.readline().strip()  # Read lines that contain actual values

            if args.verbose:
                print(f'#### Processing files: Indices of columns (name (index)): ID ({idx_snp}); REF ({idx_ref}); ALT ({idx_alt}); INFO ({idx_info})')
                print(f'#### Processing files: Loop through input file')

            count = 0
            while line != '':
                tmp_lst = line.split()
                if args.use_rsid:
                    snp = f'{tmp_lst[idx_chr]}:{tmp_lst[idx_pos]}:{tmp_lst[idx_ref]}:{tmp_lst[idx_alt]}'
                    rsid = tmp_lst[idx_rsid]
                else:
                    snp = tmp_lst[idx_snp]
                    rsid = '-'  # Assume rsID is missing
                ref = tmp_lst[idx_ref]
                alt = tmp_lst[idx_alt]
                info_val = tmp_lst[idx_info].split(';')

                # 2022/10 update: allow unordered fields
                alt_frq, maf, rsq, genotyped = -1, -1, -1, []
                for info in info_val:
                    try:
                        info_key, info_value = info.split('=')
                        if info_key == args.col_names[0]: #AF
                            alt_frq = info_value
                        elif info_key == args.col_names[1]: #MAF
                            maf = info_value
                        elif info_key == args.col_names[2]: #R2
                            rsq = info_value
                    except: # IMPUTED/TYPED/TYPED_ONLY, cannot be split by ';'
                        if args.mixed_genotype_status:
                            genotyped.append(info) # Allow mixed genotype status for one SNP
                        else:
                            genotyped = info # If all variants has single genotype status

                # Handle mixed genotype
                # ALL=all genotyped, SOME=at least one genotyped, NONE=no genotyped, -1=field missing
                if args.mixed_genotype_status:
                    genotyped = genotype_status(genotyped, args.genotyped, args.imputed)

                # Sanity check
                # Raise error if any field of AF, MAF, R2 and genotype status is missing
                for field in [alt_frq, maf, rsq, genotyped]:
                    if field == -1:
                        print("\nError: Missing field(s) in AF, MAF, Rsq, Genotyped status (missing value indicated by -1):")
                        print(f"  - SNP={snp}, rsID={rsid}, REF={ref}, ALT={alt}, AF={alt_frq}, MAF={maf}, ImputationQuality(R2)={rsq}, GenotypeStatus={genotyped}")
                        print("Exit")
                        exit()

                output_fh.write(f'{snp}\t{rsid}\t{ref}\t{alt}\t{alt_frq}\t{maf}\t{rsq}\t{genotyped}\n')
                count += 1
                if args.verbose:
                    if count % 1000 == 0:
                        print('.', end='', flush=True)
                    if count % 50000 == 0:
                        print(f'{count} lines processed')
                line = input_fh.readline().strip()
        output_fh.close()
        if args.verbose:
            print(f'\n#### Processing files: Done. input={args.input[i]}')

if __name__ == "__main__":
    write_info('')



