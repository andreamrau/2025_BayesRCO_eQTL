#!/bin/sh
#$ -l h_vmem=10G
#$ -e /travail/araul/GS_data/error
#$ -o /travail/araul/GS_data/log
#$ -q longq

# Create cluster file
awk -F' ' '{print $1" "$2" "$1}' WGS_300_GENE-SWitCH.filter_maf0.05_geno0.1.fam > cluster.txt

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
  FULLCHR=chr${CHR}
  # Filtering by chromosomes and remove missing genotypes
  plink --bfile  WGS_300_GENE-SWitCH.filter_maf0.05_geno0.1 --chr ${CHR} --geno 0 --out WGS_300_GS_${FULLCHR} --make-bed 
  # Provide unique names for variants
  awk 'BEGIN{FS=OFS="\t"} NR==0{print;next} {print $1, "SNP"NR, $3, $4, $5, $6}' WGS_300_GS_${FULLCHR}.bim > temp
  mv temp WGS_300_GS_${FULLCHR}.bim
  # Calculate per-breed allele frequency stats and HWE stats for variant filtering
  plink --bfile  WGS_300_GS_${FULLCHR} --freq --within cluster.txt --out WGS_300_GS_${FULLCHR}
  plink --bfile  WGS_300_GS_${FULLCHR} --hardy --out WGS_300_GS_${FULLCHR}
  ## Identify SNPS to keep after per-breed MAF filter (5%) and removing fixed heterozygotes
  Rscript allchr_preprocessing_geno_data_noprune_perbreed-filters.R ${FULLCHR}
  ## Filter down to retained variants
  plink --bfile WGS_300_GS_${FULLCHR} --extract WGS_300_GS_${FULLCHR}_keep.txt --out WGS_300_GS_filter_${FULLCHR} --make-bed

done




