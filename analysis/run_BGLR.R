library(BGLR)
genodir <- "/travail/araul/GS_data/"
exprdir <- "/travail/araul/GS_data/dat/"
split <- choose <- function(x, split = "_", pos = 2) {
     tmp <- unlist(lapply(strsplit(x, split = split, fixed=TRUE),
                   function(xx) xx[pos]))
     return(tmp)
}

## Loop over all expression files ---------------------------------------------------------
for(dat_file in list.files(exprdir)) {

  for(annot in c("none", "vep", "atacseq.methylation3.matched", "vep.atacseq.methylation3.matched")) {
    
  cat(dat_file, "\n")
  
  ## Analysis parameters --------------------------------------------------------------------
  learning <- split(split(dat_file, pos=2), split=".", pos=2)
  chr <- split(dat_file, pos=5)
  ensembl <- split(dat_file, pos=6)
  gene <- split(dat_file, pos=7)
  tissue <- split(split(dat_file, pos=8), split=".", pos=1)
  model <- "BGLR"
  SAVEPATH <- paste0("/travail/araul/GS_data/BGLR_results/", annot, "_results/")
  DATPATH <- "/travail/araul/GS_data/"
  BASE <- paste(chr, gene, tissue, model, learning, sep="_")

  ## Learning genotypes, recode to 0,1,2 format ---------------------------------------------
  system(paste0("cp ", DATPATH, "WGS_300_GS_filter_", chr, ".fam ", SAVEPATH, BASE, ".fam"))
  system(paste0("cp ", DATPATH, "WGS_300_GS_filter_", chr, ".bed ", SAVEPATH, BASE, ".bed"))
  system(paste0("cp ", DATPATH, "WGS_300_GS_filter_", chr, ".bim ", SAVEPATH, BASE, ".bim"))
  genotypes <- read_bed(bed_file = paste0(SAVEPATH, BASE, ".bed"),
                        bim_file = paste0(SAVEPATH, BASE, ".bim"),
                        fam_file = paste0(SAVEPATH, BASE, ".fam"))
  out <- genotypes$x
  out[out==2] <- NA
  out[out==3] <- 2
  X <- matrix(out, ncol=genotypes$n, nrow=genotypes$p, byrow=TRUE)
  X <- t(X)
  X <- scale(X, center=TRUE, scale=TRUE)

  rm(genotypes)
  rm(out)
  gc()

  ## Subset variants by annotation
  if(annot == "none") {
    X <- X
  } else {
    var_subset <- read.table(paste0("/travail/araul/GS_data/WP2_annotations/", annot, "/", annot,
                                    "_annot_WGS_GS_filter_", tissue, "-stages123_", chr, ".txt"))
    var_choice <- which(rowSums(var_subset[,-ncol(var_subset)]) > 0)
    X <- X[,var_choice]
  }

  ## Gene expression ------------------------------------------------------------------------
  expr <- read.table(paste0(exprdir, dat_file))
  fullexpr <- read.table(paste0(SAVEPATH, BASE, ".fam"), stringsAsFactors = FALSE)
  expr_ids <- read.table(paste0(DATPATH, learning, ".keep"), stringsAsFactors = FALSE)
  fullexpr$V6 <- NA
  fullexpr[match(expr_ids$V2, fullexpr$V2),"V6"] <- expr

  ## Provice GRM and use Bayesian RKHS
  G <- tcrossprod(X)
  G <- G/ncol(X)  # G <- G/mean(diag(G))
  EVD <- eigen(G)
  gblup <- BGLR(y=fullexpr$V6,
                ETA=list(list(V=EVD$vectors, d=EVD$values, model='RKHS')),
                saveAt=paste0(SAVEPATH, BASE, "_rkhs_"), nIter=12000, burnIn=2000)
  varE <- gblup$varE
  varU <- gblup$ETA[[1]]$varU
  ## Heritability
  h2 <- varU/(varU+varE)
  h2

  ## Predicted values
  gv <- data.frame(breed=fullexpr$V1, ID=fullexpr$V2, yHat=gblup$yHat)
  yTrue <- rep(NA, nrow(gv))
  names(yTrue) <- gv$ID
  DU_ids <- read.table(paste0(DATPATH, "DU.keep"), stringsAsFactors = FALSE)$V2
  LD_ids <- read.table(paste0(DATPATH, "LD.keep"), stringsAsFactors = FALSE)$V2
  LW_ids <- read.table(paste0(DATPATH, "LW.keep"), stringsAsFactors = FALSE)$V2
  DU_expr <- LD_expr <- LW_expr <- dat_file
  substr(DU_expr, 9, 10) <- "DU"
  substr(LD_expr, 9, 10) <- "LD"
  substr(LW_expr, 9, 10) <- "LW"
  yTrue[match(DU_ids, names(yTrue))] <- unlist(read.table(paste0(exprdir, DU_expr)))
  yTrue[match(LD_ids, names(yTrue))] <- unlist(read.table(paste0(exprdir, LD_expr)))
  yTrue[match(LW_ids, names(yTrue))] <- unlist(read.table(paste0(exprdir, LW_expr)))
  gv$yTrue <- yTrue

  rm(gblup)
  gc()
  
  ## Save results and clean-up
  write.table(h2, file=paste0(SAVEPATH, BASE, "_h2.txt"), row.names=FALSE,
              col.names=FALSE, quote=FALSE)
  write.table(gv, file=paste0(SAVEPATH, BASE, "_gv.txt"), row.names=FALSE,
              col.names=TRUE, sep="\t", quote=FALSE)

  system(paste0("rm ", SAVEPATH, BASE, ".fam"))
  system(paste0("rm ", SAVEPATH, BASE, ".bed"))
  system(paste0("rm ", SAVEPATH, BASE, ".bim"))
}
}
