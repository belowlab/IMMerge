for study in  MEGA_AA MEGA_Asn MEGA_HA MEGA_AI MEGA_other; do
	for chr in `seq 1 22` X; do
		echo "SNP	REF(0)	ALT(1)	Genotyped	ALT_Frq	MAF	Rsq"  > /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info
                 echo "module add samtools
		bcftools query -f'%ID\t%REF\t%ALT\tIMPUTED\t%AF\t%MAF\t%R2_HAT\n'  /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz >> /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info
                gzip -f /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info " > infoJobs/${study}.${chr}.sh
		sbatch --job-name="${chr}${study}" --mem=2g -t 03:00:00 --wrap="sh infoJobs/${study}.${chr}.sh"
	done
done

#		echo "awk 'BEGIN {OFS=\"\\t\"}; {if (NR==1) print \"SNP\tREF(0)\tALT(1)\tGenotyped\tALT_Frq\tMAF\tRsq\"; else if (\$8>1) print \$1,\$4,\$5,\"Imputed\",\$6,\$7,\"1\"; else print \$1,\$4,\$5,\"Imputed\",\$6,\$7,\$8}' /proj/epi/CVDGeneNas/hhighlan/PAGE/PAGE3/updateVCFplayground/${study}/chr${chr}.anno2 >  /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.info


# MEGA_HI didn't have anyone in topmed so the info file used above doesn't exist
for study in MEGA_HI; do
	for chr in `seq 1 22` X; do 
		echo "SNP	REF(0)	ALT(1)	Genotyped	ALT_Frq	MAF	Rsq"  > /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info
		 echo "module add samtools
		bcftools query -f'%ID\t%REF\t%ALT\tIMPUTED\t%AF\t%MAF\t%R2\n'  /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz >> /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info
		gzip -f /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info " > infoJobs/${study}.${chr}.sh
		sbatch --job-name="${chr}${study}" --mem=1g -t 03:00:00 --wrap="sh infoJobs/${study}.${chr}.sh"


		#bcftools view -H /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/MEGA_HI/MEGA_HI.chr22.vcf.gz |cut -f3,4,5,8 | tr \; \\t | cut -f 1,2,3,4,5,6,7 | sed 's/MAF=//g;s/AF=//g;s/R2=//g' | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,"0.5i",$6,$7,"-\t-\t-\t-\t-"}' >> /proj/epi/Genetic_Data_Center/calico/PAGEIII/GWAS_TOPMedFreeze8Imputed/finalized/${study}/${study}.chr${chr}.vcf.gz.info
	done
done
