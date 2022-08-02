# IMMerge

## Required packages and versions
1. This project works with python 3.7 and above. Below packages are needed:
	* python modules:
	  * pandas 1.3.3
	  * xopen 1.4.0
	* Commandline tool:
	  * bgzip (Can be downloaded from https://github.com/samtools/htslib)
2. Packages that might need manual installation: pandas, bgzip
	1. To install missing packages use: ```pip install package_name```
	2. For example: ```pip install xopen```

## Arguments
1. The program is designed to be used in command line with user input flags and values.
2. Valid flags are:
	* ```--input```: (Required) files to be merged, multiple files are allowed
	* ```--info```:  (Optional) Directory/name to info files. Default path is the same directory as corresponding input file, default info file share the same name as input file, except for suffix (.info.gz)
	* ```--output```: (Optional) Default is merged.vcf.gz and saved at current working directory. Output file name without suffix.
	* ```--thread```: (Optional) Defines how many threads to use in multiprocessing.
		* Default value is 1.
		* Valid values are integers. If number of threads <0, will use 1 instead of user supplied value
	* ```--missing```: (Optional, default is 0) Defines number of missing values allowed for each variant.
		* Cannot exceed total number of files to be merged.
		* If a variant is missing for some individuals, the values will be ".|." (or other user supplied value with --na_rep as "na_rep|na_rep") in merged output file.
		* If --missing is 0, only variants shared by all input files will be saved in merged result.
	* ```--na_rep```: (Optional) Defines what symbol to use for missing values. Default is ".". This flag is ignored if --missing is 0.
	* ```--r2_threshold```: (Optional, default is 0, ie. no filtering) Only variants with combined imputation quality score r2_combined≥r2_threshold will be saved in the merged file
	* ```--r2_output```: Default is "z_transformation". Defines how imputation quality score is calculated in the output file. Valid values are:
		* z_transformation **(recommended)**: Fisher z-transformation
		* weighted_average: calculated weighted average of r2. Weight is determined by number of individuals of each file.
		* mean, min, max: Mean, min or max of r2, ignore missing values.
		* first: output imputation quality score r2 from the first file. In order to use this setting and avoid missing r2, "--missing" must be 0.
	* ```--r2_cap```: (Optional) Adjust imputation quality score r2 by --r2_cap if imputation quality score=1. Only valid for z transformation to avoid infinity
	* ```--duplicate_id```: (Optional, default is 0) Defines number of duplicate individuals in each input file (usually as a sanity check of imputation in subsiles). Duplicate IDs should be the first N columns in each file and not mixed with unique IDs.
		* For example, the first 100 individuals in each input file are duplicated on purpose. Set --duplicated_id to 100 so that only the first set of these IDs will be saved in output file.
		* (ie. starting from the second input file, data of the first 100 individuals will be skipped in the merged output)
	* ```--check_duplicate_id```: (Optional) Default is False. Check if there are duplicate IDs, then rename non-first IDs to ID:2, ID:3, ..., ID:index_of_input_file+1.
	* ```--write_with```: (Optional) Default is bgzip. Write to bgziped file with bgzip. User can supply specific path to bgzip such as ```/user/bin/bgzip```.
	* ```--meta_info```: (Optional) Valid values are {index of input file (1-based), 'none', 'all'}. Indicates what meta information (lines start with '##') to include in output file. Default is 1 (meta information from the first input file).
	* ```--help```: (Optional) Exit program after printing out help info, ignore any other flags and values provided.

## Calculation of combined r2 and MAF
1. r2
	1. Mean: ignore missing values in calculation
	2. Weighted average, ignore missing values in calculation
		* $$r^2_{combined} = \left( \sum_{i=1}^{n}\ r_i^2 * N_i \right) / \sum_{i=1}^{n}N_i$$
		* $r^2_i$: Imputation quality r squared (Rsq) of the i-th input file
		* $N_i$: Number of individuals in the N-th input file
		* Ignore missing values. For example a variant has below r2 in each input file:
			* File #1: Rsq=0.3, number of individuals = 1000
			* File #2: Missing, number of individuals = 2000 (← Ignore this file then)
			* File #3: Rsq=0.2, number of individuals = 3000
			* Weighted Rsq = (0.3*1000 + 0.2*3000)/(1000 + 3000) = 0.225
	3. Fisher z-transformation:
		* Adjust imputation quality score as $R^2 = R^2 - r2_cap$
		* z-transformation: $z = \frac{1}{2}ln\frac{1+r}{1-r}$
		* Take weighted average of z
		* Convert z back to R using tanh function: $R = \frac{e^z - e^{-z}}{e^z + e^{-z}}$
		* Square R to get combined imputation quality
2. MAF: weighted average, ignore missing values: Use the same equation as weighted Rsq
	* $$MAF_{combined} = \left( \sum_{i=1}^{n}MAF_i * N_i \right) / \sum_{i=1}^{n}N_i $$
	* $MAF_i$: Minor allele frequency of the i-th input file
	* $N_i$: Number of individuals in the N-th input file

## File format and special notes:
1. All input files should be compressed (gzip or bgzip). File format should follow post-imputation VCF file from TOMed.
2. Variants in input vcf.gz files should be sorted by position (TOPMed output is already sorted).
3. Corresponding compressed info files should be stored in the same directory as [file_name].info.gz, unless individually specified using ```--info``` flag.
4. Use ```make_info.py``` script if .info.gz files are missing. Input VCFs must at leat have these four fields in the INFO column: AF (alt allele frequency), MAF (Minor allele frequency), imputation quality score (such as R2), Genotyped (such as IMPUTED/TYPED/TYPED_ONLY)
5. Do not move or modify variants_retained.info.txt and variants_excluded.info.txt until current run is completed.
6. Output files will be overwritten if another run saves output in the same directory with the same file name.
7. Value smaller than 0.000001 (6 digits of precision) will be rounded to 0 when outputting ALT_frq, MAF and R2 into variant_kept.txt and variant_excluded.txt. These values will also be used to replace INFO column in merged .vcf.gz file. Precision can be changed with float_format parameter in get_SNP_list.py.

## Example code
1. Example making info.gz file from input vcf.gz files. Outputs are saved in the same directory as input VCFs.
```bash
cd IMMerge
python src/merge_files.py \
	--input data_sample/sample_group1.dose.vcf.gz data_sample/sample_group2.dose.vcf.gz data_sample/sample_group3.dose.vcf.gz \
```

2. Example using sample data in ```./data_sample/```, output files are saved in ```./output_sample/```.
```bash
cd IMMerge
python src/merge_files.py \
	--input data_sample/sample_group1.dose.vcf.gz data_sample/sample_group2.dose.vcf.gz data_sample/sample_group3.dose.vcf.gz \
	--output output_sample/merged_sample \
	--check_duplicate_id true \
	--missing 1 \
	--duplicate_id 5
```
