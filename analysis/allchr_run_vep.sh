# 1 - Download the cache on genologin (where vep106 is installed)
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/release-106/variation/indexed_vep_cache/sus_scrofa_vep_106_Sscrofa11.1.tar.gz
tar xzf sus_scrofa_vep_106_Sscrofa11.1.tar.gz

# 2 - Convert plink files to vcf in interactive session, compress (and index if necessary)
# tabix -p vcf input.vcf.gz => index
srun --pty bash
module load bioinfo/PLINK/1.90b7

plink --bfile WGS_300_GS_filter_chr1 --recode vcf --out WGS_300_GS_filter_chr1
plink --bfile WGS_300_GS_filter_chr2 --recode vcf --out WGS_300_GS_filter_chr2
plink --bfile WGS_300_GS_filter_chr3 --recode vcf --out WGS_300_GS_filter_chr3
plink --bfile WGS_300_GS_filter_chr4 --recode vcf --out WGS_300_GS_filter_chr4
plink --bfile WGS_300_GS_filter_chr5 --recode vcf --out WGS_300_GS_filter_chr5
plink --bfile WGS_300_GS_filter_chr6 --recode vcf --out WGS_300_GS_filter_chr6
plink --bfile WGS_300_GS_filter_chr7 --recode vcf --out WGS_300_GS_filter_chr7
plink --bfile WGS_300_GS_filter_chr8 --recode vcf --out WGS_300_GS_filter_chr8
plink --bfile WGS_300_GS_filter_chr9 --recode vcf --out WGS_300_GS_filter_chr9
plink --bfile WGS_300_GS_filter_chr10 --recode vcf --out WGS_300_GS_filter_chr10
plink --bfile WGS_300_GS_filter_chr11 --recode vcf --out WGS_300_GS_filter_chr11
plink --bfile WGS_300_GS_filter_chr12 --recode vcf --out WGS_300_GS_filter_chr12
plink --bfile WGS_300_GS_filter_chr13 --recode vcf --out WGS_300_GS_filter_chr13
plink --bfile WGS_300_GS_filter_chr14 --recode vcf --out WGS_300_GS_filter_chr14
plink --bfile WGS_300_GS_filter_chr15 --recode vcf --out WGS_300_GS_filter_chr15
plink --bfile WGS_300_GS_filter_chr16 --recode vcf --out WGS_300_GS_filter_chr16
plink --bfile WGS_300_GS_filter_chr17 --recode vcf --out WGS_300_GS_filter_chr17
plink --bfile WGS_300_GS_filter_chr18 --recode vcf --out WGS_300_GS_filter_chr18

gzip WGS_300_GS_filter_chr1.vcf -c > WGS_300_GS_filter_chr1.vcf.gz
gzip WGS_300_GS_filter_chr2.vcf -c > WGS_300_GS_filter_chr2.vcf.gz
gzip WGS_300_GS_filter_chr3.vcf -c > WGS_300_GS_filter_chr3.vcf.gz
gzip WGS_300_GS_filter_chr4.vcf -c > WGS_300_GS_filter_chr4.vcf.gz
gzip WGS_300_GS_filter_chr5.vcf -c > WGS_300_GS_filter_chr5.vcf.gz
gzip WGS_300_GS_filter_chr6.vcf -c > WGS_300_GS_filter_chr6.vcf.gz
gzip WGS_300_GS_filter_chr7.vcf -c > WGS_300_GS_filter_chr7.vcf.gz
gzip WGS_300_GS_filter_chr8.vcf -c > WGS_300_GS_filter_chr8.vcf.gz
gzip WGS_300_GS_filter_chr9.vcf -c > WGS_300_GS_filter_chr9.vcf.gz
gzip WGS_300_GS_filter_chr10.vcf -c > WGS_300_GS_filter_chr10.vcf.gz
gzip WGS_300_GS_filter_chr11.vcf -c > WGS_300_GS_filter_chr11.vcf.gz
gzip WGS_300_GS_filter_chr12.vcf -c > WGS_300_GS_filter_chr12.vcf.gz
gzip WGS_300_GS_filter_chr13.vcf -c > WGS_300_GS_filter_chr13.vcf.gz
gzip WGS_300_GS_filter_chr14.vcf -c > WGS_300_GS_filter_chr14.vcf.gz
gzip WGS_300_GS_filter_chr15.vcf -c > WGS_300_GS_filter_chr15.vcf.gz
gzip WGS_300_GS_filter_chr16.vcf -c > WGS_300_GS_filter_chr16.vcf.gz
gzip WGS_300_GS_filter_chr17.vcf -c > WGS_300_GS_filter_chr17.vcf.gz
gzip WGS_300_GS_filter_chr18.vcf -c > WGS_300_GS_filter_chr18.vcf.gz

module purge

# 3 - module load vep in interactive session
srun --mem=20G --pty bash
module load bioinfo/Ensembl-VEP/106.1

# 4 - run
#sbatch -J vep --mem=20G --wrap="vep -i INPUT.vcf.gz -o OUTPUT.annot.merged.vcf.gz --cache --dir_cache ../vep_cache/ --species gallus_gallus --offline  --vcf --verbose --compress_output gzip --merged"
vep -i WGS_300_GS_filter_chr1.vcf.gz -o WGS_300_GS_filter_chr1_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose 
vep -i WGS_300_GS_filter_chr2.vcf.gz -o WGS_300_GS_filter_chr2_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose 
vep -i WGS_300_GS_filter_chr3.vcf.gz -o WGS_300_GS_filter_chr3_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr4.vcf.gz -o WGS_300_GS_filter_chr4_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr5.vcf.gz -o WGS_300_GS_filter_chr5_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr6.vcf.gz -o WGS_300_GS_filter_chr6_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr7.vcf.gz -o WGS_300_GS_filter_chr7_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr8.vcf.gz -o WGS_300_GS_filter_chr8_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr9.vcf.gz -o WGS_300_GS_filter_chr9_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr10.vcf.gz -o WGS_300_GS_filter_chr10_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr11.vcf.gz -o WGS_300_GS_filter_chr11_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr12.vcf.gz -o WGS_300_GS_filter_chr12_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr13.vcf.gz -o WGS_300_GS_filter_chr13_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr14.vcf.gz -o WGS_300_GS_filter_chr14_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr15.vcf.gz -o WGS_300_GS_filter_chr15_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr16.vcf.gz -o WGS_300_GS_filter_chr16_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr17.vcf.gz -o WGS_300_GS_filter_chr17_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose
vep -i WGS_300_GS_filter_chr18.vcf.gz -o WGS_300_GS_filter_chr18_vep.txt --cache --dir_cache /home/arau/.vep/ --species sus_scrofa --offline --verbose

module purge

# 5 - extract annotation column
module load statistics/R/4.3.0
R
source("allchr_clean_vep.R")
module purge










