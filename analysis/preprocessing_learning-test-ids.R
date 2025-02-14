fam <- read.table("WGS_300_GENE-SWitCH.filter_maf0.05_geno0.1.fam", header=FALSE)

LW <- fam[which(fam[,1] %in% c("LargeWhite")), 1:2]
LD <- fam[which(fam[,1] %in% c("Landrace")), 1:2]
DU <- fam[which(fam[,1] %in% c("Duroc")), 1:2]

write.table(LW, file = "LW.keep", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(LD, file = "LD.keep", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(DU, file = "DU.keep", quote=FALSE, row.names=FALSE, col.names=FALSE)



