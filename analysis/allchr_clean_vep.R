library(tidyr)
library(dplyr)
library(data.table)

## Create classes NSC and REG classes
NSC_cat <- c("coding_sequence_variant", "frameshift_variant", "inframe_deletion",
      "inframe_insertion", "initiator_codon_variant", "missense_variant",
      "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
      "stop_gained", "stop_lost", "stop_retained_variant")
REG_cat <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant", 
      "mature_miRNA_variant", "non_coding_transcript_variant", "non_coding_exon_variant", 
      "upstream_gene_variant")

for(chr_choice in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) {

  chr <- paste0("chr", chr_choice)
  cat("***", chr, "***\n")
 
  ## Read in VEP results
  vep_orig <- fread(paste0("WGS_300_GS_filter_", chr, "_vep.txt"), header = TRUE, skip = 31) %>%
    as.data.frame()

  ## Create a column for each Consequence, 1 if variant annotated and 0 otherwise
  vep <- vep_orig %>%
    select(1, Consequence) %>%
    filter(! Consequence %in% c("intron_variant", "intergenic_variant")) %>%
    unique() %>%
    separate_rows(Consequence, sep=",") %>%
    unique() %>%
    filter(Consequence %in% c(NSC_cat, REG_cat)) %>%
    mutate(class = ifelse(Consequence %in% NSC_cat, "NSC", "REG")) %>%
    mutate(annotate = 1) %>%
    select(-Consequence) %>%
    unique() %>%
    pivot_wider(names_from=class, id_cols=1, values_from=annotate) 
  colnames(vep)[1] <- "SNP" 

  ## Need to match up so same order as genotype data
  ngeno <- as.numeric(unlist(strsplit(system(paste0("wc -l WGS_300_GS_filter_", chr, ".bim"), intern=TRUE), 
    split = " "))[1])
  fullgeno <- data.frame(SNP = paste0("SNP", 1:ngeno))
  vep_final <- left_join(fullgeno, vep, by="SNP") %>%
    replace(is.na(.), 0)
  rownames(vep_final) <- vep_final[,1]
  vep_final <- vep_final[,-1]
 
  ## Write out anntotations
  write.table(vep_final, paste0("WGS_300_GS_filter_", chr, "_vep_REG-NSC.txt"),
              col.names=TRUE, row.names=FALSE, quote=FALSE)
}