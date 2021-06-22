# This code is used to merge imputed files from TopMed
# (From Hannah:)
# The goal is to combine our three post-imputation vcf.gz files
# (with around 25,000 samples in each and ~38M variants)
# into one vcf.gz file separated by chromosome
# File dir is: /data100t1/share/BioVU/TOPMed_imputation_072020 (Confirm with Hannah)

import gzip

# ----------------------- Helper functions -----------------------

# This function reads in a file with one snp in each line,
# then return a list of all SNPs.
# For this version, the first header line 'SNP' need to be skipped
def get_snp_lst(snp_fn):
    lst_snp = []
    with open(snp_fn, 'r') as fh:
        line = fh.readline()  # Skip the header line
        line = fh.readline().strip()
        while line != '':
            lst_snp.append(line)
            line = fh.readline().strip()
    return lst_snp


# Merge 3 imputed files together
def merge_files(chromosome='chr1',
                input_dir='/data100t1/share/BioVU/TOPMed_imputation_072020/',
                snp_to_keep_dir='/data100t1/home/wanying/imputation/snp_keep_or_exclude_list/',
                output_dir='/data100t1/home/wanying/imputation/merged_gz_files/'):
    print('Start merging files of '+chromosome+':')
    lst_fn = [chromosome + '_group1.dose.vcf.gz',
              chromosome + '_group2.dose.vcf.gz',
              chromosome + '_group3.dose.vcf.gz']
    snp_keep_fn = chromosome + '_merge_keep_snps.txt'

    output_fn = chromosome + '_merged.dose.vcf.gz'

    # Read in SNPs need to be kept from file
    lst_snp_to_keep = get_snp_lst(snp_to_keep_dir + snp_keep_fn)

    fh_group1 = gzip.open(input_dir + lst_fn[0], 'rt')
    fh_group2 = gzip.open(input_dir + lst_fn[1], 'rt')
    fh_group3 = gzip.open(input_dir + lst_fn[2], 'rt')
    fh_output = gzip.open(output_dir + output_fn, 'wt')  # Open for writing (appending)

    # Skip the first 18 lines (header lines of imputation info).
    # Real data starts from line 19, but need to write the header lines into merged file once
    for i in range(18):
        line_group1 = fh_group1.readline().strip()
        line_group2 = fh_group2.readline().strip()
        line_group3 = fh_group3.readline().strip()
        fh_output.write(line_group3 + '\n')

    # Need to write the header line of columns into file
    # Need to merge group1, group2 and group3 column headers
    line_group1 = fh_group1.readline().strip()
    line_group2 = fh_group2.readline().strip()
    line_group3 = fh_group3.readline().strip()
    line_to_write = line_group1 + '\t' + line_group2.split(maxsplit=109)[-1] + '\t' + line_group3.split(maxsplit=109)[
        -1]

    fh_output.write(line_to_write + '\n')

    # Read in actual data line for the while loop to start
    line_group1 = fh_group1.readline().strip()
    line_group2 = fh_group2.readline().strip()
    line_group3 = fh_group3.readline().strip()

    count = 0   # For console output
    while line_group1 != '' and line_group2 != '' and line_group3 != '':
        # when end of file is not reached for all three files
        lst_g1 = line_group1.split()
        # Remove the first 109 columns by setting maxsplit in split()
        # Those columns are repeated (with 100 duplicated individuals)
        lst_g2 = line_group2.split(maxsplit=109)
        lst_g3 = line_group3.split(maxsplit=109)
        # Get snp ID of each group and compare
        snp_group1 = lst_g1[2]
        snp_group2 = lst_g2[2]
        snp_group3 = lst_g3[2]

        # Compare snp ids with the snp_keep list file (lst_snp_to_keep)
        # Pop snps from the list until it's empty
        if lst_snp_to_keep != []:
            snp_to_keep = lst_snp_to_keep.pop(0)
            while snp_group1 != snp_to_keep:
                # To see if current sno in the vcf.dose.gz file needs to be kept
                # Continue read the file until find the snp
                line_group1 = fh_group1.readline().strip()
                lst_g1 = line_group1.split()
                snp_group1 = lst_g1[2]
            while snp_group2 != snp_to_keep:
                line_group2 = fh_group2.readline().strip()
                lst_g2 = line_group2.split(maxsplit=109)
                snp_group2 = lst_g2[2]
            while snp_group3 != snp_to_keep:
                line_group3 = fh_group3.readline().strip()
                lst_g3 = line_group3.split(maxsplit=109)
                snp_group3 = lst_g3[2]

            # Write into output file
            # Piece together 3 files (first 100 individuals are replicated)
            line_to_write = line_group1 + '\t' + lst_g2[-1] + '\t' + lst_g3[-1]

            # print('group1:\t\tline len = ', len(lst_g1))
            # print('group2:\t\tline len = ', len(lst_g2[-1].split()))
            # print('group3:\t\tline len = ', len(lst_g3[-1].split()))
            # print(len(lst_g3))
            # quit()
            # print('line to write:\tlen =',len(line_to_write.split()), '\n-------------')

            fh_output.write(line_to_write + '\n')

        # Below is to keep console busy
        count = count + 1
        if count % 100 == 0:
            print('.', end='', flush=True)

        if count % 10000 == 0:
            print('\n' + chromosome + ':', count, 'SNPs have been merged', flush=True)
        elif count % 1000 == 0:
            print('', flush=True)

        line_group1 = fh_group1.readline().strip()
        line_group2 = fh_group2.readline().strip()
        line_group3 = fh_group3.readline().strip()

    print('\n' + chromosome + ':', count, 'SNPs have been merged')
    print('\n------------------------------------\nDone')

# ----------------------- End of helper functions -----------------------
#
# chr_number_list = ['1','2','3','4','5','6','7','13','14','15','16','17','18','19','20','22','X']
#
# for i in chr_number_list:
#     check_imputation_parameters(chromosome='chr'+i)



merge_files(chromosome='chr5')
