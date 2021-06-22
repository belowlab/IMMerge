# This code is to check if imputation was performed with the correct settings
# For now this code only checks if r2 filter = 0.1
# Only need to run once

import gzip
import sys
from process_args import process_args

# File name can be passed to this code in terminal, or use import this code as in a script (need to modify a little)
args = sys.argv
verbose=True
# test_flag = True # !!! Remove this when done
dict_flags = process_args(args) # Process terminal input
if dict_flags['--verbose']!='true': verbose=False

# if test_flag:
#     print(dict_flags)

# ----------------------- Helper functions -----------------------
#
# # This function check imputation setting info based on header lines of .dose.vcf.gz files
# def check_imputatilson_parameters(chromosome='chr1', input_dir=''):
#     lst_fn = [chromosome + '_group1.dose.vcf.gz',
#               chromosome + '_group2.dose.vcf.gz',
#               chromosome + '_group3.dose.vcf.gz']
#     try:
#         fh_group1 = gzip.open(input_dir + lst_fn[0], 'rt')
#         fh_group2 = gzip.open(input_dir + lst_fn[1], 'rt')
#         fh_group3 = gzip.open(input_dir + lst_fn[2], 'rt')
#
#         # Skip the first 6 lines
#         # r2 filter data starts at line 7
#         for i in range(6):
#             fh_group1.readline()
#             fh_group2.readline()
#             fh_group3.readline()
#         # Get r2 setting
#         r2_group1 = fh_group1.readline().strip().split('=')[1]
#         r2_group2 = fh_group2.readline().strip().split('=')[1]
#         r2_group3 = fh_group3.readline().strip().split('=')[1]
#
#         if r2_group1 == r2_group2 and r2_group3 == r2_group2:
#             print(chromosome, 'r2 filter =', r2_group1)
#         else:
#             print('r2 filter is not the same across all groups!')
#             print(' - Group1 r2 filter =', r2_group1)
#             print(' - Group2 r2 filter =', r2_group2)
#             print(' - Group3 r2 filter =', r2_group3)
#
#         # close file handles
#         fh_group1.close()
#         fh_group2.close()
#         fh_group3.close()
#     except:
#         print(chromosome, 'Error: File(s) not found')
#
#
# # ----------------------- End of helper functions -----------------------
#
# chr_number_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
#                    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']
#
# for i in chr_number_list:
#     check_imputation_parameters(chromosome='chr' + i)
