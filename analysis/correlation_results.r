split.choose <- function(x, split = "_", pos = 1) {
    unlist(lapply(strsplit(x, split = split, fixed=TRUE), function(y) y[pos]))
}

datpath <- paste0("/faige_save/GS_data/dat/")
respath <- paste0("/faige_save/GS_data/acrossbreed_results/")

all.results <- vector("list", 4)
names(all.results) <- c("none", "vep", "atacseq.methylation3.matched",
                        "vep.atacseq.methylation3.matched")

true.files <- list.files(datpath)[grep(".dat$", list.files(datpath))]

for(an in names(all.results)) {

        cat("***", an, "***\n")

        datanpath <- paste0(respath, an, "_results/")
        res.files <- list.files(datanpath)
        if(length(grep(".zip$", res.files))) {
            res.files <- res.files[-c(grep(".zip$", res.files))]
        }
        if(length(grep("fort.21", res.files))) {
            res.files <- res.files[-c(grep("fort.21", res.files))]
        }

        gv.files <- res.files[grep(".gv", res.files)]
        validation.files <- gv.files[grep("validation", gv.files)]
        
        all.chr <- split.choose(validation.files, split = "_", pos = 1)
        all.genes <- split.choose(validation.files, split = "_", pos = 2)
        all.tissue <- split.choose(validation.files, split = "_", pos = 3)
        all.learning <- split.choose(validation.files, split = "_", pos = 6)
        all.validation <- split.choose(split.choose(validation.files, split = "_",
                                                    pos = 8), split=".", pos = 1)

        results <- data.frame(Chromosome=all.chr, EnsemblID = rep(NA, length(all.chr)),
                           Gene.Name=all.genes, tissue=all.tissue, annotation=an,
                           learning = all.learning,
                           validation = all.validation,
                           cor_pearson = NA, cor_spearman = NA,
                           stringsAsFactors = FALSE)

        for(i in 1:nrow(results)) {
            gene <- results$Gene.Name[i]
            chr <- results$Chromosome[i]
            tissue <- results$tissue[i]
            validation <- results$validation[i]
            learning <- results$learning[i]
            ensembl <- unique(split.choose(true.files[grep(gene, true.files)],
                                           split="_", pos=6))

            true <-
                read.table(paste0(datpath,
                                  true.files[grep(paste0(validation,  "_GS_filter_",
                                                         chr, "_", ensembl, "_",
                                                         gene, "_", tissue, ".dat"), true.files)]))

            if(learning != validation) {
            est <- read.table(paste0(datanpath, chr, "_", gene, "_", tissue, "_BayesRCO_", an, "_", learning,
                                     "_validation_", validation, ".gv"))
            } else {
            est <- read.table(paste0(datanpath, chr, "_", gene, "_", tissue, "_BayesRCO_", an, "_", learning,
                                     ".gv"))
            }

            results$cor_pearson[i] <- ifelse(!is.null(est), round(cor(true, est, method="pearson"),3), NA)
            results$cor_spearman[i] <- ifelse(!is.null(est), round(cor(true, est, method="spearman"),3), NA)
            results$EnsemblID[i] <- ensembl
        }
        all.results[[an]] <- results
    }

all.results.df <- do.call("rbind", all.results)

write.table(all.results.df, paste0("/faige_save/GS_data/correlation_results_acrossbreed.txt"),
                row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE)


## BGLR results -----------------------------------------------------------------------------------
for(annot in c("none", "vep", "atacseq.methylation3.matched", "vep.atacseq.methylation3.matched")) {
  respath <- paste0("/faige_save/GS_data/BGLR_results/", annot, "_results/")
  h2.files <- list.files(respath)[grep("_h2.txt", list.files(respath))]
  gv.files <- list.files(respath)[grep("_gv.txt", list.files(respath))]

  gblup_h2 <- matrix(NA, nrow = length(h2.files), ncol = 5)
  colnames(gblup_h2) <- c("chr", "gene", "tissue", "training", "h2")
  for(i in 1:length(h2.files)) {
    f <- h2.files[i]
    h2 <- as.numeric(read.table(paste0(respath, f)))
    gblup_h2[i,] <- c(split.choose(f, pos=1), split.choose(f, pos=2), split.choose(f, pos=3),
           split.choose(f, pos=5), h2)
  }
  write.table(gblup_h2, paste0("/faige_save/GS_data/BGLR_h2_", annot, ".txt"),
              col.names=TRUE, row.names=FALSE, quote=FALSE)

  gblup_gv <- matrix(NA, nrow = length(gv.files)*3, ncol = 7)
  colnames(gblup_gv) <- c("chr", "gene", "tissue", "training", "validation", "pearson", "spearman")
  index <- 1
  for(i in 1:length(gv.files)) {
    f <- gv.files[i]
    gv <- read.table(paste0(respath, f), header=TRUE)
    chr <- split.choose(f, pos=1)
    gene <- split.choose(f, pos=2)
    tissue <- split.choose(f, pos=3)
    training <- split.choose(f, pos=5)
    for(j in c("LargeWhite", "Duroc", "Landrace")) {
      validation <- ifelse(j == "LargeWhite", "LW", ifelse(j == "Duroc", "DU", "LD"))
      gv_subset <- gv[which(gv$breed == j),]
      pearson <- cor(gv_subset$yHat, gv_subset$yTrue, method="pearson")
      spearman <- cor(gv_subset$yHat, gv_subset$yTrue, method="spearman")
      gblup_gv[index,] <- c(chr, gene, tissue, training, validation, pearson, spearman)
      index <- index + 1
    }
  }
  write.table(gblup_gv, paste0("/faige_save/GS_data/BGLR_gv_cor_", annot, ".txt"),
                               col.names=TRUE, row.names=FALSE, quote=FALSE)
}
