The analysis scripts in this repo are organized as described below.

## Preprocessing RNA-seq data

-   `allchr_select_gene.R` : Correct logCPM expression values for sex effects

## Preprocessing WGS data

-   `preprocessing_learning-test-ids.R` : Identify animals for each of the three breeds
-   `allchr_preprocessing_geno_data_noprune.sh` : Filtering by chromosomes and remove missing genotypes, Provide unique names for variants, calculate per-breed allele frequency stats and HWE stats for variant filtering, identify SNPS to keep after per-breed MAF filter (5%) and removing fixed heterozygotes, filter down to retained variants
-   `allchr_preprocessing_geno_data_noprune_perbreed-filters.R` : Identify variants to filter out based on global heterozygosity or per-breed MAF \< 5%

## Annotations

-   `allchr_run_vep.sh` : Run VEP command line tool
-   `allchr_clean_vep.R` : Format VEP annotations (REG, NSC) for use with BayesRCO
-   `allchr_preprocessing_create-annotations.R` : Format epigenetic annotations (OCR, UMR, LMR) for use with BayesRCO

## Exploratory analysis

-   `allchr_pca_plink.sh` : Perform PCA on genotype data after LD pruning

## Run QTL mapping and prediction models

-   `run_bayesRCO.sh` : fit BayesRCpi model using BayesRCO
-   `run_BGLR.R` : fit Bayesian RKHS model using BGLR

## Results

-   `colocalization.R` : Annotate BayesRCO results with PigQTLdb QTLs
-   `correlation_results.r` : Calculate correlations between predicted and estimated expression values
