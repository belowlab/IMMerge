# This code is to test each individual script
# For debugging purpose
import get_SNP_list
import os
os.chdir('/Users/wanying/Documents/Below_lab/TOPMed_imputation_ongoing/TOPMed_merge/')

# Hard code input parameters
dict_flags = dict()
dict_flags['--input'] = ['./sample_data/subset_short_chr21_group1.dose.vcf.gz',
                         './sample_data/subset_short_chr21_group2.dose.vcf.gz',
                         './sample_data/subset_short_chr21_group3.dose.vcf.gz']
dict_flags['--missing'] = 1
dict_flags['--r2_threshold'] = 0.2
dict_flags['--r2_output'] = 'weighted_average'
dict_flags['--duplicate_id'] = 10
dict_flags['--output'] = 'merged'
dict_flags['--verbose'] = 'true'
dict_flags['--na_rep'] = 'NA'


get_SNP_list.get_snp_list(dict_flags)