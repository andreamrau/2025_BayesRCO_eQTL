library(GenomicRanges)

## Read in original annotations ------------------------------------------------------------------------------
setwd("/travail/araul/GS_data/allchr_annotations/WP2_annotations/")

## ATAC-seq
atacseq <- data.frame(read.table("all.pig.atacseq_mRp.clN_peaks.gappedPeak.txt"))
colnames(atacseq) <- c("chr", "start", "end", "annotation")

## Methylation (UMR/LMR) in muscle and liver
methylation <- data.frame(read.table("SS_liver_muscle_meth_class.txt",
                                     header=TRUE, sep="\t"))
colnames(methylation) <- c("chr", "start", "end", "tissue", "stage", "class")
methylation$tissue <- ifelse(methylation$tissue == "Liver", "liver", "muscle")
methylation$annotation <- paste0(methylation$tissue, "_", methylation$stage)
methylation$chr <- substr(methylation$chr,4,10)
methylation$annotation_full <- paste0(methylation$tissue, "_",
                                      methylation$stage, "_", methylation$class)

all_annotations <- list(
    atacseq=atacseq,
    methylation3=methylation[,c("chr", "start", "end", "annotation_full")]
)

## Find overlaps ------------------------------------------------------------------------------
for(chr_choice in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18")) {

  ## Add VEP annotations to the overall list (these are the same order as the .bim files by construction)
  all_annotations$vep <- read.table(paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/WGS_300_GS_filter_", chr_choice,
                                           "_vep_REG-NSC.txt"), header=TRUE)
  
  bim <- read.table(paste0("/travail/araul/GS_data/allchr_annotations/WGS_300_GS_filter_", chr_choice, ".bim"))
  colnames(bim)[c(1,4)] <- c("chr", "pos")
  bim_GR <- GRanges(seqnames=bim$chr, ranges=IRanges(bim$pos, bim$pos))

  for(tissue_choice in c("liver", "muscle")) {

    matched <- tissue_choice
    mismatched <- ifelse(tissue_choice == "liver", "muscle", "liver")

    overlap_annotations <- list()
    # Loop over atac/methyl/vep
    for(j in names(all_annotations)) {
      if(j == "vep") {
        overlap_annotations[[j]] <- all_annotations$vep
        next;
      }
      
      annot <- all_annotations[[j]]
      # Choose appropriate tissue
      annot_tissue <- annot[grep(tissue_choice, annot$annotation),]
      annot_GR <- GRanges(seqnames = annot_tissue$chr,
                          ranges=IRanges(annot_tissue$start, annot_tissue$end))
      full_annotations <- vector("list", length = length(unique(annot_tissue$annotation)))
      names(full_annotations) <- unique(annot_tissue$annotation)
      cat(names(full_annotations), "\n")
      
      # Loop over total number of categories (x stage, x lowly/weakly methylated)
      for(i in names(full_annotations)) {
        cat("***", i, "***\n")
        annot_choose <- annot_GR[which(annot_tissue$annotation == i)]
        overlap <- data.frame(findOverlaps(bim_GR, annot_choose))$queryHits
        tmp <- rep(0, length(bim_GR))
        tmp[overlap] <- 1
        full_annotations[[i]] <- tmp
      }
      overlap_annotations[[j]] <- do.call("cbind", full_annotations)
    }

    # none (these are constant across tissues)
    none <- rep(1, nrow(bim))
    write.table(none,
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/none/none_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    # vep (these are constant across tissues)
    vep <- overlap_annotations$vep
    vep <- data.frame(vep, other=ifelse(rowSums(vep) == 0, 1, 0))
    write.table(vep,
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/vep/vep_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(colnames(vep),
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/vep/vep_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, "_annotation-order.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    # atacseq.methylation3.matched
    atacseq.methylation3.matched <-
      data.frame(overlap_annotations$atacseq, overlap_annotations$methylation3)
    atacseq.methylation3.matched <-
      data.frame(atacseq.methylation3.matched, other=ifelse(rowSums(atacseq.methylation3.matched) == 0, 1, 0))
    write.table(atacseq.methylation3.matched,
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/atacseq.methylation3.matched/",
               "atacseq.methylation3.matched_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(colnames(atacseq.methylation3.matched),
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/atacseq.methylation3.matched/",
               "atacseq.methylation3.matched_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, "_annotation-order.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    # vep.atacseq.methylation3.matched
    vep.atacseq.methylation3.matched <-
      data.frame(overlap_annotations$vep, overlap_annotations$atacseq, overlap_annotations$methylation3)
    vep.atacseq.methylation3.matched <-
      data.frame(vep.atacseq.methylation3.matched, other=ifelse(rowSums(vep.atacseq.methylation3.matched) == 0, 1, 0))
    write.table(vep.atacseq.methylation3.matched,
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/vep.atacseq.methylation3.matched/",
               "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(colnames(vep.atacseq.methylation3.matched),
        paste0("/travail/araul/GS_data/allchr_annotations/WP2_annotations/vep.atacseq.methylation3.matched/",
               "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
               matched, "-stages123_", chr_choice, "_annotation-order.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)    
      cat(chr_choice, tissue_choice, "\n")
  }
}



