args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

strat <- read.table(paste0("WGS_300_GS_", chr, ".frq.strat"), header=TRUE, stringsAsFactors=FALSE)
strat_wide <- reshape(data = strat, direction="wide", idvar="SNP", timevar="CLST")
strat_wide$minMAF <- apply(strat_wide[,c("MAF.Duroc", "MAF.Landrace", "MAF.LargeWhite")], 1, min)

hwe <- read.table(paste0("WGS_300_GS_", chr, ".hwe"), header=TRUE, stringsAsFactors = FALSE)
hwe <- hwe[which(hwe$TEST == "ALL"),]
het <- as.numeric(unlist(lapply(strsplit(hwe$GENO, split = "/", fixed=TRUE), function(x) x[2]) ))


remove <- data.frame(SNP = strat_wide$SNP, minMAF = ifelse(strat_wide$minMAF < 0.05, 1, 0),
                     allhet = ifelse(het == 300, 1, 0), stringsAsFactors = FALSE)
keep_index <- remove[which(remove$minMAF + remove$allhet == 0), "SNP"]
keep_index <- as.character(keep_index)

write.table(keep_index, paste0("WGS_300_GS_", chr, "_keep.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
