#!/bin/sh
#$ -l h_vmem=10G
#$ -e /travail/araul/GS_data/error
#$ -o /travail/araul/GS_data/log
#$ -q longq

cd /travail/araul/GS_data/pca/

## Only keep chromosomes 1-18
plink-1-09 --bfile WGS_300_GENE-SWitCH.filter_maf0.05_geno0.1 --chr 1-18 --out WGS_300_chr1.18 --make-bed
# Provide unique names for variants
awk 'BEGIN{FS=OFS="\t"} NR==0{print;next} {print $1, "SNP"NR, $3, $4, $5, $6}' WGS_300_chr1.18.bim > temp
mv temp WGS_300_chr1.18.bim
# Perform linkage pruning prior to PCA (identify prune sites)
plink-1-09 --bfile WGS_300_chr1.18 --indep-pairwise 50 10 0.1 --out WGS_300_chr1.18_prune
# Perform PCA 
plink-1-09 --bfile WGS_300_chr1.18 --extract WGS_300_chr1.18_prune.prune.in --make-bed --pca --out WGS_300_GS_chr1.18_prune_allchr
