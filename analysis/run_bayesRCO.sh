#!/bin/sh
#$ -l h_vmem=5G
#$ -e /travail/araul/GS_data/error
#$ -o /travail/araul/GS_data/log
#$ -q longq
#$ -m e
#$ -M andrea.rau@inrae.fr
#$ -N chr7_SUPT3H_muscle_LW_vep.atacseq.methylation3.matched
#$ -l h='!(node130|node120|node190)'
GENE="SUPT3H"
ENSEMBL="ENSSSCG00000001709"
CHR="chr7"
LEARNING="LW"
MODEL="BayesRCO"
TISSUE="muscle"
ANNOTCHOICE="vep.atacseq.methylation3.matched"
NCAT=12

## CHANGE THIS ACCORDING TO ANNOTATIONS
SAVEPATH="/travail/araul/GS_data/acrossbreed_results/${ANNOTCHOICE}_results"  
ANNOTPATH="/travail/araul/GS_data/WP2_annotations"
SOFTPATH="/travail/araul/BayesRCO/src"
DATPATH="/travail/araul/GS_data"
ANNOTATIONS="${ANNOTPATH}/${ANNOTCHOICE}/${ANNOTCHOICE}_annot_WGS_GS_filter_${TISSUE}-stages123_${CHR}.txt"

## Temporary copy of genome files
BASE=${CHR}_${GENE}_${TISSUE}_${MODEL}_${ANNOTCHOICE}_${LEARNING}

cp ${DATPATH}/WGS_300_GS_filter_${CHR}.fam ${SAVEPATH}/${BASE}.fam
cp ${DATPATH}/WGS_300_GS_filter_${CHR}.bed ${SAVEPATH}/${BASE}.bed
cp ${DATPATH}/WGS_300_GS_filter_${CHR}.bim ${SAVEPATH}/${BASE}.bim

## Keep track of a common allele to set as A1 (use full population of n=300 for this)
awk '{print $2,$5}' ${SAVEPATH}/${BASE}.bim > ${SAVEPATH}/${BASE}_A1allele.txt

## Subset to include correct learning breed and two remaining validation breeds
# Modify appropriate .fam files for learning data (validation data = NA)
/bao/plink/plink-1.09-x86_64/plink --bfile ${SAVEPATH}/${BASE} --keep ${DATPATH}/${LEARNING}.keep --reference-allele ${SAVEPATH}/${BASE}_A1allele.txt --make-bed --out ${SAVEPATH}/${BASE}_training 
awk 'FNR==NR{a[NR]=$1;next}{$6=a[FNR]}1' ${DATPATH}/dat/WGS_100.${LEARNING}_GS_filter_${CHR}_${ENSEMBL}_${GENE}_${TISSUE}.dat ${SAVEPATH}/${BASE}_training.fam > temp_${BASE}
mv temp_${BASE} ${SAVEPATH}/${BASE}_training.fam

for VALBREED in DU LD LW
  do
  /bao/plink/plink-1.09-x86_64/plink --bfile ${SAVEPATH}/${BASE} --keep ${DATPATH}/${VALBREED}.keep --reference-allele ${SAVEPATH}/${BASE}_A1allele.txt --make-bed --out ${SAVEPATH}/${BASE}_validation_${VALBREED} 
  awk '$6="NA"' ${SAVEPATH}/${BASE}_validation_${VALBREED}.fam > temp_${BASE}
  mv temp_${BASE} ${SAVEPATH}/${BASE}_validation_${VALBREED}.fam
done 

## Move to results directory
cd ${SAVEPATH}

# BayesRCO learning 
${SOFTPATH}/bayesRCO -bfile ${BASE}_training -out ${BASE} -catfile ${ANNOTATIONS} -ncat ${NCAT} \
    -ndist 5 -gpin 0.0,0.00001,0.0001,0.001,0.01 
#    -numit 20000 -burnin 0 -thin 1 
#    -ndist 6 -gpin 0.0,0.000001,0.00001,0.0001,0.001,0.01
#    -ndist 5 -gpin 0.0,0.00001,0.0001,0.001,0.01
#    -ndist 4 -gpin 0.0,0.0001,0.001,0.01

# BayesRCO validation 
for VALBREED in DU LD LW
  do
  ${SOFTPATH}/bayesRCO -bfile ${BASE}_validation_${VALBREED} -predict -out ${BASE}_validation_${VALBREED} \
  -model ${BASE}.model -freq ${BASE}.frq -param ${BASE}.param -ncat ${NCAT} -catfile ${ANNOTATIONS} \
  -ndist 5 -gpin 0.0,0.00001,0.0001,0.001,0.01 
done

## Remove temporary copy of genome files
rm ${SAVEPATH}/${BASE}.fam
rm ${SAVEPATH}/${BASE}.bed
rm ${SAVEPATH}/${BASE}.bim
rm ${SAVEPATH}/${BASE}_training.fam
rm ${SAVEPATH}/${BASE}_training.bed
rm ${SAVEPATH}/${BASE}_training.bim
rm ${SAVEPATH}/${BASE}_validation_*.fam
rm ${SAVEPATH}/${BASE}_validation_*.bed
rm ${SAVEPATH}/${BASE}_validation_*.bim
rm ${SAVEPATH}/${BASE}_A1allele.txt






