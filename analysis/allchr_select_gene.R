# library(data.table)

setwd("/travail/araul/GS_data/Pig_expression")

## Sex correction for all genes prior to global PCA
expr.muscle <- read.table("LogCPM_counts_GENE-SWitCH_muscle.txt", header=TRUE)
texpr.muscle <- t(expr.muscle)
colnames(texpr.muscle) <- unlist(strsplit(texpr.muscle[1,], split = "_M"))
texpr.muscle <- texpr.muscle[-1,]
texpr <- apply(texpr.muscle,2,as.numeric)
rownames(texpr) <- rownames(texpr.muscle)
gexpr.muscle <- texpr

expr.liver <- read.table("LogCPM_counts_GENE-SWitCH_liver.txt", header=TRUE)
texpr.liver <- t(expr.liver)
colnames(texpr.liver) <- unlist(strsplit(texpr.liver[1,], split = "_M"))
texpr.liver <- texpr.liver[-1,]
texpr <- apply(texpr.liver,2,as.numeric)
rownames(texpr) <- rownames(texpr.liver)
gexpr.liver <- texpr

meta <- read.table("/travail/araul/GS_data/pca/WGS_300_GENE-SWitCH.filter_maf0.05_geno0.1.fam")
colnames(meta)[c(1,2,5)] <- c("breed", "name", "sex")
meta$breed <- ifelse(meta$breed == "Duroc", "DU", ifelse(meta$breed == "LargeWhite", "LW", "LD"))
meta <- meta[,-c(3,4,6)]
meta$sex <- as.character(meta$sex)
meta$name <- as.character(meta$name)
meta$breed <- as.character(meta$breed)

## Muscle
cat("Muscle!\n")
pca.muscle <- matrix(NA, nrow=nrow(gexpr.muscle), ncol=300)
colnames(pca.muscle) <- colnames(gexpr.muscle)
rownames(pca.muscle) <- rownames(gexpr.muscle)
 for(i in 1:nrow(gexpr.muscle)) {
   ## Format expression data
    expr.dat <- gexpr.muscle
    gexpr <- expr.dat[i,]
    ## Correct for sex effect and multiply expression residuals by 100
    for(breed in c("DU", "LD", "LW")) {
      index <- grep(breed, colnames(expr.dat))
      if(breed == "LW") {
        expression_subset <- as.numeric(gexpr)[index]
        meta_subset <- meta[index,]
        mod <- lm(expression_subset ~ 1)
        expression_residuals <- residuals(mod)
        expression_choose <- expression_residuals * 100
        pca.muscle[i,index] <- expression_choose
      } else {
        expression_subset <- as.numeric(gexpr)[index]
        meta_subset <- meta[index,]
        mod <- lm(expression_subset ~ meta_subset$sex)
        expression_residuals <- residuals(mod)
        expression_choose <- expression_residuals * 100
        pca.muscle[i,index] <- expression_choose
      }
    }
}
write.table(pca.muscle, paste0("/travail/araul/GS_data/Pig_expression/allchr_muscle_sexcorrected.txt"),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

## Liver
cat("Liver!\n")
pca.liver <- matrix(NA, nrow=nrow(gexpr.liver), ncol=300)
colnames(pca.liver) <- colnames(gexpr.liver)
rownames(pca.liver) <- rownames(gexpr.liver)
 for(i in 1:nrow(gexpr.liver)) {
   cat(i, "\n")
   ## Format expression data
    expr.dat <- gexpr.liver
    gexpr <- expr.dat[i,]
    ## Correct for sex effect and multiply expression residuals by 100
    for(breed in c("DU", "LD", "LW")) {
      index <- grep(breed, colnames(expr.dat))
      if(breed == "LW") {
        expression_subset <- as.numeric(gexpr)[index]
        meta_subset <- meta[index,]
        mod <- lm(expression_subset ~ 1)
        expression_residuals <- residuals(mod)
        expression_choose <- expression_residuals * 100
        pca.liver[i,index] <- expression_choose
      } else {
        expression_subset <- as.numeric(gexpr)[index]
        meta_subset <- meta[index,]
        mod <- lm(expression_subset ~ meta_subset$sex)
        expression_residuals <- residuals(mod)
        expression_choose <- expression_residuals * 100
        pca.liver[i,index] <- expression_choose
      }
    }
}
write.table(pca.liver, paste0("/travail/araul/GS_data/Pig_expression/allchr_liver_sexcorrected.txt"),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

   







