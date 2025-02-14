library(dplyr)
library(tidyr)
library(ggplot2)
library(mixOmics)
library(cowplot)
library(forcats)
library(ComplexHeatmap)
library(viridis)
library(ChIPpeakAnno)
library(cowplot)
library(ggdendro)
library(GenomicRanges)
library(fgsea)
library(gghighlight)
library(rtracklayer)
library(lemon)
library(BGLR)
library(karyoploteR)
library(BSgenome)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.UCSC.susScr11.refGene)
library(ggplotify)
library(data.table)

setwd("figures/")

## Auxillary functions ---------------------------------------------------------

split_choose <- function(x, split = "_", pos = 2) {
    tmp <- unlist(lapply(strsplit(x, split = "_", fixed=TRUE),
                         function(xx) xx[pos]))
    return(tmp)
}


join_overlap_full <- function(x,y){
    options(stringsAsFactors = F)
    hit <- findOverlaps(x,y)
    colnames(mcols(x)) <- paste0(colnames(mcols(x)),".x")
    colnames(mcols(y)) <- paste0(colnames(mcols(y)),".y")
    cns.x <- mcols(x) %>% colnames()
    cns.y <- mcols(y) %>% colnames()

    z <- x[queryHits(hit),]
    hitz <- cbind(mcols(y)[subjectHits(hit),],
                  mcols(x)[queryHits(hit),])
    colnames(hitz) <- c(cns.y,cns.x)
    mcols(z) <- hitz
    uniqex <- x[-queryHits(hit),]
    uniqey <- y[-subjectHits(hit),]
    mcols(uniqex)[, (length(cns.x)+1):(length(cns.x)+length(cns.y))] <- NA
    colnames(mcols(uniqex)) <- c(cns.x,cns.y)
    mcols(uniqex) <- mcols(uniqex)[,c((length(cns.x)+1):(length(cns.x)+length(cns.y)),
                                      1:(length(cns.x)))]
    mcols(uniqey)[, (length(cns.y)+1):(length(cns.x)+length(cns.y))] <- NA
    colnames(mcols(uniqey)) <- c(cns.y,cns.x)
    res <- c(z, uniqex,uniqey)
    mcols(res) <- mcols(res) %>% data.frame() %>% mutate_all(as.character())
    return(res)
}

hyp_to_annot_enrichment = function(hyp_file, nclass = 4, annotations_name){
    nannot = length(annotations_name)
    df_enrichment = data.frame(Annotation = rep(annotations_name,each=nclass),
                               Class = rep(1:nclass, nannot),
                               Freq = rep(0, nclass*nannot))
    for (i in 1:nannot){
        df_enrichment$Freq[which(df_enrichment$Annotation == 
                                   annotations_name[i])] = colMeans(hyp_file[,paste0("Nk",1:nclass,"_",i,sep="")]/rowSums(hyp_file[,paste0("Nk",1:nclass,"_",i,sep="")]))
    }
    return(df_enrichment)
}

## Figure 2 (breed heterogeneity) ----------------------------------------------
### A: Genomic PCA -------------------------------------------------------------
eigenvec <- read.table("data/WGS_300_GS_chr1.18_prune_allchr.eigenvec")
eigenval <- read.table("data/WGS_300_GS_chr1.18_prune_allchr.eigenval")
colnames(eigenvec)[c(1:4)] <- c("breed", "ID", "PCA1", "PCA2")
perc1 <-
    data.frame(PC1 = paste0("PC1: ",round(eigenval$V1[1] / sum(eigenval$V1)*100,1), "%"))
perc2 <-
    data.frame(PC2 = paste0("PC2: ",round(eigenval$V1[2] / sum(eigenval$V1)*100,1), "%"))
perc <- cbind(perc1, perc2)
perc$perc <- paste0(perc$PC1, "\n", perc$PC2)

fig1a_tmp <- ggplot(eigenvec) +
    geom_point(aes(x=PCA1, y=PCA2, col=breed)) +
    geom_hline(aes(yintercept = 0), lty=2, col="grey30") +
    geom_vline(aes(xintercept = 0), lty=2, col="grey30") +
    geom_text(data=perc, aes(x=Inf, y=-Inf, label=perc), hjust = 1, vjust = -.1) +
    theme_bw() +
    xlab("PC1") + ylab("PC2") +
    scale_color_discrete(name="") +
    theme(legend.justification = "top")
fig1a <- fig1a_tmp + theme(legend.position="none")

legend <- get_legend(
    # create some space to the left of the legend
    fig1a_tmp + theme(legend.box.margin = margin(0, 0, 0, 8))
)

### B: Transcriptomic PCA ------------------------------------------------------
muscle_noncorrected <- read.table("data/LogCPM_counts_GENE-SWitCH_muscle.txt",
                                  header=TRUE)
liver_noncorrected <- read.table("data/LogCPM_counts_GENE-SWitCH_liver.txt",
                                 header=TRUE)
rownames(muscle_noncorrected) <-
    unlist(strsplit(muscle_noncorrected$SampleID, split="_M", fixed=TRUE))
muscle_noncorrected <- muscle_noncorrected[,-1]
rownames(liver_noncorrected) <-
    unlist(strsplit(liver_noncorrected$SampleID, split="_L", fixed=TRUE))
liver_noncorrected <- liver_noncorrected[,-1]
muscle_pca <- pca(muscle_noncorrected, scale=TRUE, ncomp=10)
liver_pca <- pca(liver_noncorrected, scale=TRUE, ncomp=10)

muscle_perc <- data.frame(tissue="muscle", PC1=muscle_pca$prop_expl_var$X[1],
                          PC2=muscle_pca$prop_expl_var$X[2])
liver_perc <- data.frame(tissue="liver", PC1=liver_pca$prop_expl_var$X[1],
                         PC2=liver_pca$prop_expl_var$X[2])
perc <- bind_rows(muscle_perc, liver_perc)
perc$perc <- paste0("PC1: ", round(perc$PC1*100,1), "%\nPC2 :", 
                    round(perc$PC2*100, 1), "%")

muscle_pca <- muscle_pca$variates$X
liver_pca <- liver_pca$variates$X
transcriptome_pca <- rbind(
    data.frame(ID = rownames(muscle_pca), tissue = "muscle", muscle_pca),
    data.frame(ID = unlist(strsplit(rownames(liver_pca), split="_L", fixed=TRUE)),
               tissue = "liver", liver_pca)) %>%
    mutate(breed = substr(ID, 1, 2))

fig1b <- ggplot(transcriptome_pca) +
    geom_point(aes(x=PC1, y=PC2, col=breed)) +
    geom_hline(aes(yintercept = 0), lty=2, col="grey30") +
    geom_vline(aes(xintercept = 0), lty=2, col="grey30") +
    geom_text(data=perc, aes(x=Inf, y=-Inf, label=perc), hjust = 1, vjust = -.1) +
    facet_wrap(~tissue) +
    theme_bw() +
    scale_color_discrete(name="") +
    theme(legend.position="none")



### C Boxplots of logCP expression --------------------------------------------

logcpm_liver <-
    t(read.table("data/LogCPM_counts_GENE-SWitCH_liver_11genes.txt")) %>%
    as.data.frame()
logcpm_muscle <-
    t(read.table("data/LogCPM_counts_GENE-SWitCH_muscle_11genes.txt")) %>%
    as.data.frame()
logcpm_liver$pop <- substr(rownames(logcpm_liver), 1, 2)
logcpm_liver$tissue <- "liver"
logcpm_muscle$pop <- substr(rownames(logcpm_muscle), 1, 2)
logcpm_muscle$tissue <- "muscle"
logcpm_liver_long <- logcpm_liver %>% as.data.frame() %>%
    gather(key = "gene", value="expression", -pop, -tissue)
logcpm_muscle_long <- logcpm_muscle %>% as.data.frame() %>%
    gather(key = "gene", value="expression", -pop, -tissue)

all.dat <- bind_rows(logcpm_liver_long, logcpm_muscle_long)

all.dat$pop <- factor(all.dat$pop)
levels(all.dat$pop) <- c("Duroc", "Landrace", "LargeWhite")
## Remove L3HYPDH because no longer included in Crespo et al
fig1c <- ggplot(all.dat %>% filter(gene != "L3HYPDH")) +
    geom_boxplot(aes(x=gene, y=expression, fill = pop, color= pop), alpha = 0.3) +
    geom_hline(aes(yintercept=0), lty = 2, col="black") +
    facet_wrap(~tissue, ncol=1) +
    theme_bw() +
    ylab("Log CPM expression") + xlab("") +
    theme(legend.position="none")


### D: Heritability density plots ----------------------------------------------
gblup_results <- vector("list", 1)
names(gblup_results) <- c("none")
gblup_h2 <- vector("list", 1)
names(gblup_h2) <- c("none")
for(i in names(gblup_results)) {
    gv <- read.table(paste0("data/BGLR_gv_cor_", i, ".txt"), header=TRUE)
    gv$annotation <- i
    gblup_results[[i]] <- gv
    h2 <- read.table(paste0("data/BGLR_h2_", i, ".txt"), header=TRUE)
    h2$annotation <- i
    gblup_h2[[i]] <- h2
}

gblup_results <- bind_rows(gblup_results)
gblup_h2 <- bind_rows(gblup_h2)
gblup_results$annotation <- factor(gblup_results$annotation)
gblup_results$annotation <- fct_relevel(gblup_results$annotation,
                                        "none")
gblup_h2$annotation <- factor(gblup_h2$annotation)
gblup_h2$annotation <- fct_relevel(gblup_h2$annotation,
                                   "none")
levels(gblup_h2$annotation) <- c("none")

## Remove L3HYPDH because no longer included in Crespo et al
fig1d <- ggplot(gblup_h2 %>% filter(annotation == "none", gene != "L3HYPDH")) +
    geom_density(aes(x=h2, fill=training, color=training), alpha = 0.3) +
    geom_rug(aes(x=h2, color=training), size=1.5, alpha=0.5) +
    facet_wrap(~tissue, ncol=2) +
    xlab(expression(paste("h"^"2"))) +
    theme_bw() +
    theme(legend.position="none")

### E: Learning validations  ---------------------------------------------------

tmp <- full_join(gblup_results %>% filter(training == validation),
                 gblup_h2, by= c("chr", "gene", "tissue", "training", "annotation"))
## Remove L3HYPDH because no longer included in Crespo et al
fig1e <- ggplot(tmp %>% filter(annotation == "none", gene != "L3HYPDH")) +
    geom_boxplot(aes(x=training, y=spearman, col=training, fill=training), alpha =0.5) +
    geom_jitter(aes(x=training, y=spearman, col=training), alpha =0.2, width=0.1) +
    xlab("Training breed") +
    ylab("Correlation") +
    guides(col="none", fill="none") +
    facet_wrap(~tissue) +
    theme_bw()

### Put it all together --------------------------------------------------------
fig1top <-
    plot_grid(fig1a, fig1b, legend, rel_widths=c(1,2, .4), labels=c("A","B", ""),
              ncol=3)
fig1mid <- plot_grid(fig1c, labels=("C"))
fig1bot <-
    plot_grid(fig1d, fig1e, rel_widths=c(1,1), labels=c("D","E"))
fig2 <- plot_grid(fig1top, fig1mid, fig1bot,
                  rel_heights = c(1,2,1), ncol=1)

## Figure 3 (annotations) ------------------------------------------------------

## Read in data
ensembl <- rtracklayer::import('data/susScr11.ensGene.gtf.gz')
refGene <- rtracklayer::import('data/susScr11.refGene.gtf.gz')
df_ensembl = as.data.frame(ensembl)
df_ensembl$strand = "*"
gr_ensembl = makeGRangesFromDataFrame(df_ensembl, keep.extra.columns = T)
df_cat <- read.table("data/Traits_category.txt",
                     header=T,
                     sep="\t")

annot_bim <- read.table("data/allchr_annotations.bim", header=T)
annot_file <- read.table("data/allchr_annotations.txt",header=T)
annot_file$Annotation = rep("Not annotated", nrow(annot_file))
annot_file$Annotation[which(annot_file$REG == 1 | annot_file$NSC ==1)] = "REG+NSC"
annot_file$Annotation[which(annot_file$REG == 0 & annot_file$NSC == 0 &
                                annot_file$other == 0 &
                                annot_file$tissue == "liver")]= "liver-specific"
annot_file$Annotation[which(annot_file$REG == 0 & annot_file$NSC == 0 &
                                annot_file$other == 0 &
                                annot_file$tissue == "muscle")]= "muscle-specific"
QTLdb <- read.table("data/QTLdb_pigSS11.bed",
                    skip = 18,
                    fill = T,
                    row.names = NULL,
                    header = F,
                    sep="\t")
QTLdb = QTLdb[,-c(5:12)]
colnames(QTLdb) = c("chromosome","start","end","trait")
QTLdb$trait = sub("\\(.*", "", QTLdb$trait)
QTLdb$chromosome = sub("....","", QTLdb$chromosome)
QTLdb$full_trait = gsub(" QTL ", "",QTLdb$trait)
QTLdb = left_join(QTLdb,
                  df_cat)
gr_QTLdb = makeGRangesFromDataFrame(QTLdb, keep.extra.columns = T)
gr_QTLdb_filtered = makeGRangesFromDataFrame(filter(as.data.frame(gr_QTLdb), width <10),
                                             keep.extra.columns = T)

### A: pie charts ------------------------------------------------------------------

df_freq_annot = 
  data.frame(SNP = rep(c("Annotated", "Non annotated"),2),
             Frequency = c(table(filter(annot_file, tissue=="muscle")$other)/
                             nrow(filter(annot_file, tissue=="muscle")),
                                         table(filter(annot_file, tissue=="liver")$other)/
                             nrow(filter(annot_file, tissue=="liver"))),
             Tissue = rep(c("Muscle + REG + NSC","Liver + REG + NSC"),each=2))

df_freq_annot$pos <-
    c(0.92, .4, 0.9, 0.4)
df_freq_annot$Tissue <- factor(df_freq_annot$Tissue)
levels(df_freq_annot$Tissue) <- c("Liver", "Muscle")
df_freq_annot$SNP <- factor(df_freq_annot$SNP)
levels(df_freq_annot$SNP) <- c("VEP and epigenetic annotations", "Unannotated")
fig3a <- ggplot(df_freq_annot, aes(x="", y=Frequency, fill=SNP)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+
    geom_text(aes(y=pos, x=1.2, label=paste0(100*round(Frequency, 3), "%")),
              col="white", fontface = 2) +
    facet_rep_wrap(~Tissue)+
    theme_void()+ xlab("") + ylab("") +
    scale_fill_manual(values = c("black", "grey60"), name = "") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.position="bottom", strip.text.x = element_text(size = 12))


### B: overlaps ------------------------------------------------------------------

df_SNP_size = rbind(cbind(as.data.frame(table(rowSums(filter(annot_file, 
                                                             tissue=="muscle")[3:13]))),
                          Set="Muscle + REG + NSC"),
                    cbind(as.data.frame(table(rowSums(filter(annot_file, 
                                                             tissue=="liver")[3:13]))),
                          Set="Liver + REG + NSC"))
colnames(df_SNP_size)[1:2]=c("nSNPs","Size")

options(scipen=10000)
df_SNP_size_plot <- filter(df_SNP_size, nSNPs!=0)
df_SNP_size_plot$Set <- factor(df_SNP_size_plot$Set)
levels(df_SNP_size_plot$Set) <- c("liver",
                                  "muscle")
fig3b <- ggplot(df_SNP_size_plot, aes(x=nSNPs, y=Size, fill=Set))+
    geom_bar(stat="identity", position=position_dodge2())+
    # scale_y_continuous(trans='log10')+
    theme_bw()+
    theme( strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank())+
    xlab("Number of annotations")+
    ylab("Number of SNPs")+
    scale_y_sqrt() +
    scale_fill_manual(name="", values = c("#414487FF", "#22A884FF")) +
    theme(legend.position="inside", legend.position.inside = c(.65,.65))

options(scipen=10)

### C: per annotation -------------------------------------------------------------------

df_annot_size = data.frame(Annotation = rep(colnames(annot_file)[3:13],2),
                           Size = c(colSums(filter(annot_file, 
                                                   tissue=="muscle")[3:13]),
                                    colSums(filter(annot_file, 
                                                   tissue=="liver")[3:13])),
                           Tissue = rep(c("muscle","liver"),each=11))

df_annot_size$Tissue[grep("REG",df_annot_size$Annotation)] = "not tissue-specific"
df_annot_size$Tissue[grep("NSC",df_annot_size$Annotation)] = "not tissue-specific"

df_annot_size$Annotation2 = rep(c("REG", "NSC",
                                  "OCR 30dpf", "OCR 70dpf", "OCR newborn",
                                  "UMR 30dpf", "LMR 30dpf",
                                  "UMR 70dpf", "LMR 70dpf",
                                  "UMR newborn", "LMR newborn"),2)
df_annot_size$Annotation2 = factor(df_annot_size$Annotation2,
                                   levels= c("REG", "NSC",
                                             "OCR 30dpf", "OCR 70dpf", "OCR newborn",
                                             "UMR 30dpf", "UMR 70dpf", "UMR newborn",
                                             "LMR 30dpf", "LMR 70dpf", "LMR newborn"))

fig3c <- ggplot(df_annot_size, aes(x= Annotation2, y= Size, fill=Tissue))+
    geom_bar(stat="identity", position =position_dodge())+
    coord_flip() +
    theme_bw()+
    theme(
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=10),
          axis.title = element_text(size=12),
          strip.text = element_text(size=12), legend.position="bottom",
          legend.justification = .8)+
    scale_fill_manual(name = "", values = c("#414487FF", "#22A884FF", "black")) +
    guides(color = "none") +
    # guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    xlab("") + ylab("Number of SNPs")


### D: dendrogram ------------------------------------------------------------------

temp = cbind(filter(annot_file, tissue == "muscle")[5:13],
             filter(annot_file, tissue == "liver")[5:13])
colnames(temp)=c(paste("muscle_",
                       colnames(annot_file)[5:13],
                       sep=""),
                 paste("liver_",
                       colnames(annot_file)[5:13],
                       sep=""))
temp = cbind(temp, filter(annot_file, tissue == "muscle")[3:4])

dist_jaccard = dist(t(temp), method="binary")

hclust_jaccard = hclust(dist_jaccard)

cl <- cutree(hclust_jaccard,k=20)

dhc <- as.dendrogram(hclust_jaccard)

labels_x <- names(cl[hclust_jaccard$order])
dat <- dendro_data(hclust_jaccard)$labels

ddata <- dendro_data(dhc, type = "rectangle")
label_data <- bind_cols(filter(segment(ddata), x == xend & x%%1 == 0 & yend ==0),
                        label=names(cl[hclust_jaccard$order]))
label_data$Tissue = "Not tissue specific"
label_data$Tissue[grep("liver",label_data$label)]="Liver"
label_data$Tissue[grep("muscle",label_data$label)]="Muscle"
label_data$label2 = c("REG",
                      "LMR 30dpf","LMR 70dpf","LMR newborn",
                      "LMR 30dpf","LMR 70dpf","LMR newborn",
                      "NSC",
                      "UMR 30dpf","UMR 70dpf",
                      "UMR 70dpf","UMR newborn",
                      "UMR 30dpf", "UMR newborn",
                      "OCR newborn","OCR 30dpf","OCR 70dpf",
                      "OCR newborn","OCR 30dpf","OCR 70dpf")

# Rectangular lines
fig3d <- ggplot(segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() +
    scale_y_reverse(expand = c(1, 0))+
    theme_dendro()+
    geom_point(data = dat, aes(x = x, y = y,color = label_data$Tissue)) +
    geom_text(data=label_data, aes(x=xend, y=yend-0.01, label=label2, hjust=0,
                                   color = label_data$Tissue), size=4)+
    guides(color=guide_legend(title="Tissue"))+
    scale_color_manual(name = "", values = c("#414487FF", "#22A884FF", "black")) +
    theme(legend.position="none" ,
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

### E: Fisher ------------------------------------------------------------------

annot_pos = data.frame(seqnames = annot_bim$chr,
                       start = annot_bim$bp_pos,
                       end = annot_bim$bp_pos,
                       strand = NA,
                       allele_1 = annot_bim$allele_1,
                       allele_2 =  annot_bim$allele_2)


df_annot_muscle = cbind(annot_pos,
                        filter(annot_file, tissue=="muscle"))
df_annot_liver = cbind(annot_pos,
                       filter(annot_file, tissue=="liver"))

gr_annot_muscle = makeGRangesFromDataFrame(df_annot_muscle,
                                           keep.extra.columns = T,
                                           seqnames.field="seqnames")
gr_annot_liver = makeGRangesFromDataFrame(df_annot_liver,
                                          keep.extra.columns = T,
                                          seqnames.field="seqnames")

gr_full_muscle = join_overlap_full(gr_annot_muscle,gr_QTLdb_filtered)
gr_full_filtered_muscle = subsetByOverlaps(subsetByOverlaps(gr_full_muscle,
                                                            gr_annot_muscle), 
                                           gr_QTLdb_filtered)

gr_left_muscle = subsetByOverlaps(gr_full_muscle, gr_annot_muscle)
df_left_muscle =as.data.frame(gr_left_muscle)
df_left_muscle$cat_trait.y[is.na(df_left_muscle$cat_trait.y)]="None"
df_left_muscle$QTL = "Yes"
df_left_muscle$QTL[which(df_left_muscle$cat_trait.y == "None")]="No"

gr_full_liver = join_overlap_full(gr_annot_liver,gr_QTLdb_filtered)
gr_full_filtered_liver = subsetByOverlaps(subsetByOverlaps(gr_full_liver,
                                                           gr_annot_liver), 
                                          gr_QTLdb_filtered)

gr_left_liver = subsetByOverlaps(gr_full_liver, gr_annot_liver)
df_left_liver =as.data.frame(gr_left_liver)
df_left_liver$cat_trait.y[is.na(df_left_liver$cat_trait.y)]="None"
df_left_liver$QTL = "Yes"
df_left_liver$QTL[which(df_left_muscle$cat_trait.y == "None")]="No"

list_fisher_reg = sapply(13:14,
                         function(x) {fisher.test(table(df_left_muscle$QTL,
                                                        df_left_muscle[,x]))})

df_fisher_reg = data.frame(pval = apply(list_fisher_reg, 2, 
                                        function(x) {x$p.value}),
                           OR = apply(list_fisher_reg, 2, 
                                      function(x) {x$estimate}),
                           IC_min = apply(list_fisher_reg, 2, 
                                          function(x) {x$conf.int[1]}),
                           IC_max = apply(list_fisher_reg, 2, 
                                          function(x) {x$conf.int[2]}),
                           Tissue = rep("",2),
                           Modality = c("REG","NSC"),
                           Stage = c("",""))

list_fisher_muscle = sapply(15:23,
                            function(x) {fisher.test(table(df_left_muscle$QTL,
                                                           df_left_muscle[,x]))})

df_fisher_muscle = data.frame(pval = apply(list_fisher_muscle, 2, 
                                           function(x) {x$p.value}),
                              OR = apply(list_fisher_muscle, 2, 
                                         function(x) {x$estimate}),
                              IC_min = apply(list_fisher_muscle, 2, 
                                             function(x) {x$conf.int[1]}),
                              IC_max = apply(list_fisher_muscle, 2, 
                                             function(x) {x$conf.int[2]}),
                              Tissue = rep("Muscle",9),
                              Modality = c("OCR","OCR","OCR",
                                           "UMR","LMR",
                                           "UMR","LMR",
                                           "UMR","LMR"),
                              Stage = c("30dpf", "70dpf", "newborn",
                                        "30dpf","30dpf",
                                        "70dpf","70dpf",
                                        "newborn","newborn"))

list_fisher_liver = sapply(15:23,
                           function(x) {fisher.test(table(df_left_liver$QTL,
                                                          df_left_liver[,x]))})

df_fisher_liver = data.frame(pval = apply(list_fisher_liver, 2, 
                                          function(x) {x$p.value}),
                             OR = apply(list_fisher_liver, 2, 
                                        function(x) {x$estimate}),
                             IC_min = apply(list_fisher_liver, 2, 
                                            function(x) {x$conf.int[1]}),
                             IC_max = apply(list_fisher_liver, 2, 
                                            function(x) {x$conf.int[2]}),
                             Tissue = rep("Liver",9),
                             Modality = c("OCR","OCR","OCR",
                                          "UMR","LMR",
                                          "UMR","LMR",
                                          "UMR","LMR"),
                             Stage = c("30dpf", "70dpf", "newborn",
                                       "30dpf","30dpf",
                                       "70dpf","70dpf",
                                       "newborn","newborn"))

df_fisher_all = rbind(df_fisher_reg,
                      df_fisher_muscle,
                      df_fisher_liver)

df_fisher_all$Modality = factor(df_fisher_all$Modality,
                                levels=c("LMR","UMR","OCR","NSC","REG"))
df_fisher_all$Stage = factor(df_fisher_all$Stage,
                             levels=c("","newborn","70dpf","30dpf"))
fig3e <- ggplot(df_fisher_all,
       aes(y=interaction(Modality,Stage, sep = " "), x=OR, label=Stage, 
           color = -log10(pval))) +
    geom_point(size=4) +
    geom_errorbarh(aes(xmin=IC_min, xmax=IC_max), height=.3) +
    geom_vline(xintercept=1, linetype='longdash') +
    facet_wrap(~Tissue, ncol=1, scale="free_y")+
    scale_color_gradient(low = "#FCBBA1", high="#A50F15", name = "-log10 p-value")+
    xlab("Odds Ratio") + ylab("") +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position="bottom",
          legend.justification = .8)


### Put it all together --------------------------------------------------------
fig3left <-
    plot_grid(fig3a, fig3b, fig3c, fig3d,
              labels=c("A","B", "C", "D"), ncol=2)
# 958 x 737
fig3 <- plot_grid(fig3left, fig3e, labels = c("", "E"), rel_widths = c(2,1))


## Figure 4 (QTL enrichment) ---------------------------------------------------

df_cat <- read.table("data/Traits_category.txt",
                     header=T,
                     sep="\t")

fgsea_res_bayesR <- read.table("data/df_bayesR_list_fgsea_scenarios.csv",
                               sep=",",
                               header=T)
fgsea_res_bayesR=fgsea_res_bayesR[,-1]
fgsea_res_bayesR$gene = gsub("_.*", "\\1", fgsea_res_bayesR$scenario)
fgsea_res_bayesR$tissue = "liver"
fgsea_res_bayesR$tissue[grep("muscle", fgsea_res_bayesR$scenario)] = "muscle"
fgsea_res_bayesR$learning = gsub(".*_", "\\1", fgsea_res_bayesR$scenario)
fgsea_res_bayesR$QTL_trait = gsub(" QTL ", "",fgsea_res_bayesR$pathway)
fgsea_res_bayesR = left_join(fgsea_res_bayesR,
                             df_cat,
                             by=c("QTL_trait"="full_trait"))

fgsea_res_bayesRCO <- read.table("data/df_bayesRCO_list_fgsea_scenarios.csv",
                                 sep=",",
                                 header=T)
fgsea_res_bayesRCO=fgsea_res_bayesRCO[,-1]
fgsea_res_bayesRCO$gene = gsub("_.*", "\\1", fgsea_res_bayesRCO$scenario)
fgsea_res_bayesRCO$tissue = "liver"
fgsea_res_bayesRCO$tissue[grep("muscle", fgsea_res_bayesRCO$scenario)] = "muscle"
fgsea_res_bayesRCO$learning = gsub(".*_", "\\1", fgsea_res_bayesRCO$scenario)
fgsea_res_bayesRCO$QTL_trait = gsub(" QTL ", "",fgsea_res_bayesRCO$pathway)
fgsea_res_bayesRCO = left_join(fgsea_res_bayesRCO,
                               df_cat,
                               by=c("QTL_trait"="full_trait"))

fgsea_res_models = rbind(cbind(fgsea_res_bayesR, model = "bayesR"),
                         cbind(fgsea_res_bayesRCO, model= "bayesRCO"))

fgsea_res_models_cat = read.table("data/df_bayes_models_list_fgsea_scenarios_categories.csv",
                                  sep=",",
                                  header=T)
fgsea_res_models_cat$Freq = rep(1, nrow(fgsea_res_models_cat))
fgsea_res_models_cat_0.05 = filter(fgsea_res_models_cat, pval <0.05)

### A: model enrichment ------------------------------------------------------------------

## Remove L3HYPDH because no longer included in Crespo et al
fig4a <- ggplot(fgsea_res_models %>% filter(gene != "L3HYPDH"), 
                aes(x=-log10(pval), fill= model, color=model))+
    geom_density(alpha=0.5)+
    theme_bw()+
    scale_fill_manual(values=c("grey30", "#F8765CFF"), name="") +
    scale_color_manual(values=c("grey30", "#F8765CFF"), name="") +
    facet_wrap(~tissue) +
    xlab("-log10 p-value") +
    theme(legend.position="none", legend.position.inside = c(.8,.8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

## KS liver
x <- fgsea_res_models %>% 
  filter(tissue == "liver", model == "bayesR", gene != "L3HYPDH") %>% 
  dplyr::select(pval)
y <- fgsea_res_models %>% 
  filter(tissue == "liver", model == "bayesRCO", gene != "L3HYPDH") %>% 
  dplyr::select(pval)
ks.test(unlist(x),unlist(y))

## KS muscle
x <- fgsea_res_models %>% 
  filter(tissue == "muscle", model == "bayesR", gene != "L3HYPDH") %>% 
  dplyr::select(pval)
y <- fgsea_res_models %>% 
  filter(tissue == "muscle", model == "bayesRCO", gene != "L3HYPDH") %>% 
  dplyr::select(pval)
ks.test(unlist(x),unlist(y))

### 4B: category enrichment ------------------------------------------------------------------

## Remove L3HYPDH because no longer included in Crespo et al
fgsea_res_models_cat_0.05 = filter(fgsea_res_models_cat, pval <0.05, gene != "L3HYPDH")

df_freq_tissue_cat = aggregate(fgsea_res_models_cat_0.05$Freq,
                               list(fgsea_res_models_cat_0.05$model,
                                    fgsea_res_models_cat_0.05$tissue,
                                    fgsea_res_models_cat_0.05$pathway),
                               sum)

colnames(df_freq_tissue_cat)=c("Model","Tissue","Trait category","Freq")
df_freq_tissue_cat$`Trait category` <- factor(df_freq_tissue_cat$`Trait category`)
levels(df_freq_tissue_cat$`Trait category`)[3] <- "Meat and carcass"
tmp <- data.frame(Model = c("BayesRCO", "BayesR", "BayesR"),
                  Tissue = c("liver", "muscle", "muscle"),
                  `Trait category` = c("Exterior", "Production",
                                       "Meat and carcass"),
                  Freq = c(NA, NA, NA), check.names=FALSE)
tmp2 <- bind_rows(df_freq_tissue_cat, tmp)
tmp2$Model <- factor(tmp2$Model)
levels(tmp2$Model) <- c("BayesR", paste0("BayesRC", "\u03C0"))
fig4b <- ggplot(tmp2,
                aes(x=`Trait category`, y=Freq,
                                        fill=Model, color=Model))+
    geom_bar(stat="identity", position=position_dodge(preserve = "single"),
             alpha = 0.5)+
    facet_grid(~Tissue)+
    scale_fill_manual(values=c("grey30", "#F8765CFF"), name="") +
    scale_color_manual(values=c("grey30", "#F8765CFF"), name="") +
    theme_bw()+ xlab("") + ylab("Frequency") +
    coord_flip()


### C: gene enrichment ---------------------------------------------------------
fgsea_res_bayesRCO$cat_trait <- factor(fgsea_res_bayesRCO$cat_trait)
levels(fgsea_res_bayesRCO$cat_trait)[3] <- "Meat and carcass"
## Remove L3HYPDH because no longer included in Crespo et al
fig4c <- ggplot(fgsea_res_bayesRCO %>% filter(gene != "L3HYPDH"), 
                aes(y=-log10(pval),x=gene))+
    facet_rep_grid(tissue~cat_trait)+
    geom_boxplot(alpha=0.3)+
    theme_bw()+
    ylab("-log10 p-value") +
    xlab("") +
    theme(strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=10),
          axis.title = element_text(size=14),
          strip.text = element_text(size=12))+
    coord_flip()+
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2)

### Put it all together --------------------------------------------------------
fig4top <-
    plot_grid(fig4a, fig4b, labels=c("A","B"), ncol=2)
# 880 x 729
fig4 <- plot_grid(fig4top, fig4c, labels = c("", "C"),
                  rel_heights = c(1,2), ncol=1)

## Figure 5 (across-breed predictions) ------------------------------------------

results <- read.table("data/correlation_results_acrossbreed.txt", header=TRUE)
results$annotation <- factor(results$annotation)
results$annotation <-
    fct_relevel(results$annotation, "none", "vep",
                "atacseq.methylation3.matched","vep.atacseq.methylation3.matched")
levels(results$annotation) <- c("none", "vep", "atac+methyl", "vep+atac+methyl")

### A: heatmap ------------------------------------------------------------------
## Remove L3HYPDH because no longer included in Crespo et al
bayesrco_results <- results %>%
    filter(annotation == "vep+atac+methyl", Gene.Name != "L3HYPDH")
bayesrco_summary <- bayesrco_results %>%
    dplyr::select(-cor_pearson, -Chromosome, -EnsemblID, -Gene.Name, -annotation) %>%
    mutate(cor_spearman = ifelse(learning == validation, NA, cor_spearman)) %>%
    group_by(learning, validation, tissue) %>%
    summarize(average_spearman = mean(cor_spearman),
              SD_spearman = sd(cor_spearman)) %>%
    ungroup()
bayesrco_summary_learning <- bayesrco_results %>%
    dplyr::select(-cor_pearson, -Chromosome, -EnsemblID, -Gene.Name, -annotation) %>%
    filter(learning != validation) %>%
    group_by(learning, validation, tissue) %>%
    summarize(average_spearman = mean(cor_spearman)) %>%
    ungroup()

fig5a <- ggplot(bayesrco_summary) +
    geom_tile(aes(x=learning, y=validation,, fill=average_spearman)) +
    geom_text(data = bayesrco_summary_learning,
              aes(x=learning, y=validation, label=round(average_spearman, 4))) +
    facet_wrap(~tissue, ncol=1) +
    scale_fill_continuous(type = "gradient", low = "white", high="#B63679FF",
                          na.value="grey90", name="Correlation") +
    theme_minimal() +
    xlab("Learning breed") + ylab("Validation breed") +
    theme(legend.position = "right")

### B: correlation -------------------------------------------------------------

## Remove L3HYPDH because no longer included in Crespo et al
results_perf <- results %>%
    mutate(combo = paste0(learning, "\u2192", validation)) %>%
    filter(learning != validation,
           annotation %in% c("vep+atac+methyl"), Gene.Name != "L3HYPDH")
results_perf$combo <- factor(results_perf$combo)
fig5b <- ggplot(results_perf) +
    geom_point(aes(x=Gene.Name, y=cor_spearman, color=combo,
                   shape=combo), size = 3, alpha = 0.6) +
    facet_wrap(~tissue, nrow = 2, strip.position="right") +
    geom_hline(yintercept = 0, lty = 2) +
    theme_bw() +
    ylim(c(-.3,1)) +
    xlab("") + ylab("Validation correlation") +
    scale_color_manual(name="",
                         values = c("#f8766d", "#f8766d",
                                    "#00BA38", "#00BA38",
                                    "#619CFF", "#619CFF")) +
    scale_shape_manual(name="", values = c(19, 17, 17, 19, 17, 19)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    coord_flip()

### C: breed enrichment --------------------------------------------------------
df_cat <- read.table("data/Traits_category.txt",
                     header=T,
                     sep="\t")
fgsea_res_bayesRCO <- read.table("data/df_bayesRCO_list_fgsea_scenarios.csv",
                                 sep=",",
                                 header=T)
fgsea_res_bayesRCO=fgsea_res_bayesRCO[,-1]
fgsea_res_bayesRCO$gene = gsub("_.*", "\\1", fgsea_res_bayesRCO$scenario)
fgsea_res_bayesRCO$tissue = "liver"
fgsea_res_bayesRCO$tissue[grep("muscle", fgsea_res_bayesRCO$scenario)] = "muscle"
fgsea_res_bayesRCO$learning = gsub(".*_", "\\1", fgsea_res_bayesRCO$scenario)
fgsea_res_bayesRCO$QTL_trait = gsub(" QTL ", "",fgsea_res_bayesRCO$pathway)
fgsea_res_bayesRCO = left_join(fgsea_res_bayesRCO,
                               df_cat,
                               by=c("QTL_trait"="full_trait"))
fgsea_res_bayesRCO$cat_trait <- factor(fgsea_res_bayesRCO$cat_trait)
levels(fgsea_res_bayesRCO$cat_trait)[3] <- "Meat and carcass"

## Remove L3HYPDH because no longer included in Crespo et al
fig5c <- ggplot(fgsea_res_bayesRCO %>% filter(gene != "L3HYPDH"), 
                aes(y=-log10(pval),x=cat_trait, fill=learning))+
    geom_boxplot(alpha=0.6)+
    theme_bw()+ ylab("-log10 p-value") + xlab("") +
    scale_fill_discrete(name="Learning breed")+
    facet_wrap(~tissue, ncol=1, strip.position="right") +
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
    coord_flip() +
    theme(legend.position="bottom")

### D: window v ------------------------------------------------------------------

Vbeta_annot <- read.table("data/Vbeta_annot_all.txt",
                          header=T,
                          sep="\t")
Vbeta_annot = Vbeta_annot[-grep("tissue",Vbeta_annot$tissue),]
Vbeta_annot$start = as.numeric(Vbeta_annot$start)
Vbeta_annot$end = as.numeric(Vbeta_annot$end)
Vbeta_annot$Sum_vi = as.numeric(Vbeta_annot$Sum_vi)

df_cor_Vbeta = data.frame(Scenario = rep(NA, 66),
                          Cor_Vbeta = rep(NA, 66),
                          Cor_Vbeta_sum = rep(NA, 66))
compt = 1
for (i in unique(Vbeta_annot$gene)){
    for (j in unique(Vbeta_annot$tissue)){
        temp = filter(Vbeta_annot, gene== i & tissue==j)
        df_cor_Vbeta$Scenario[compt] = paste0(i, "_", j, "_DU_LD",sep="")
        df_cor_Vbeta$Scenario[compt+1] = paste0(i, "_", j, "_DU_LW",sep="")
        df_cor_Vbeta$Scenario[compt+2] = paste0(i, "_", j, "_LW_LD",sep="")
        df_cor_Vbeta$Cor_Vbeta_sum[compt] =
          cor(as.numeric(filter(temp, learning=="DU")$Sum_vi),
              as.numeric(filter(temp, learning=="LD")$Sum_vi))
        df_cor_Vbeta$Cor_Vbeta_sum[compt+1] = 
          cor(as.numeric(filter(temp, learning=="DU")$Sum_vi),
              as.numeric(filter(temp, learning=="LW")$Sum_vi))
        df_cor_Vbeta$Cor_Vbeta_sum[compt+2] = 
          cor(as.numeric(filter(temp, learning=="LW")$Sum_vi),
              as.numeric(filter(temp, learning=="LD")$Sum_vi))
        df_cor_Vbeta$Cor_Vbeta[compt] = 
          cor(as.numeric(filter(temp, learning=="DU")$Vbeta),
              as.numeric(filter(temp, learning=="LD")$Vbeta))
        df_cor_Vbeta$Cor_Vbeta[compt+1] = 
          cor(as.numeric(filter(temp, learning=="DU")$Vbeta),
              as.numeric(filter(temp, learning=="LW")$Vbeta))
        df_cor_Vbeta$Cor_Vbeta[compt+2] = 
          cor(as.numeric(filter(temp, learning=="LW")$Vbeta),
              as.numeric(filter(temp, learning=="LD")$Vbeta))
        compt = compt + 3
    }
}

df_cor_Vbeta$Gene = rep(unique(Vbeta_annot$gene), each = 6)
df_cor_Vbeta$Tissue = rep(unique(Vbeta_annot$tissue), each = 3)
df_cor_Vbeta$Breed1 = rep(c("DU","DU","LW"))
df_cor_Vbeta$Breed2 = rep(c("LD","LW","LD"))
df_cor_Vbeta$Pairwise = rep(c("DU_LD","DU_LW","LW_LD"))
df_cor_Vbeta$Scenario1 = paste0(df_cor_Vbeta$Gene,
                                "_",
                                df_cor_Vbeta$Tissue,
                                "_",
                                df_cor_Vbeta$Breed1,
                                sep="")
df_cor_Vbeta$Scenario2 = paste0(df_cor_Vbeta$Gene,
                                "_",
                                df_cor_Vbeta$Tissue,
                                "_",
                                df_cor_Vbeta$Breed2,
                                sep="")
df_cor_Vbeta$Pairwise <- factor(df_cor_Vbeta$Pairwise)
levels(df_cor_Vbeta$Pairwise) <- c("(DU,LD)", "(DU,LW)", "(LW,LD)")

## Remove L3HYPDH because no longer included in Crespo et al
fig5d <- ggplot(df_cor_Vbeta %>% filter(Gene != "L3HYPDH"),
                aes(x=Cor_Vbeta_sum, color = Pairwise,fill=Pairwise))+
    geom_density(alpha=0.3)+
    theme_bw()+
    facet_wrap(~Tissue)+
    xlab(expression(paste("Window-based ", V[i], " correlation")))+
    xlim(c(.8,1))+
    scale_colour_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"),
                        name="")+
    scale_fill_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"),
                      name="") +
    theme(legend.position="bottom")

### Put it all together --------------------------------------------------------
fig5top <-
    plot_grid(fig5a, fig5b, rel_widths=c(1,1.75),
              labels=c("A","B"), ncol=2)
fig5bot <-
    plot_grid(fig5c, fig5d, rel_widths=c(1,1),
              labels=c("C","D"), ncol=2)
fig5 <- plot_grid(fig5top, fig5bot, ncol=1, rel_heights=c(1,1))


## Figure 6 (IGF2) -------------------------------------------------------------

gene_choice <- "IGF2"
ensembl_choice <- "ENSSSCG00000035293"
chr_choice <- "chr2"
tissue_choice <- "liver"
window <- 1000000 
pos_choice <- 50750725

### A0: karyotype -----------------------------------------------

txdb <- TxDb.Sscrofa.UCSC.susScr11.refGene
all.genes <- genes(txdb)
head(all.genes)
gene_labels <- data.frame(chr="chr2", pos = 1469132,
                          label = "IGF2")
region <- GRanges(data.frame(seqnames = "chr2",
                     start = pos_choice-window,
                     end = pos_choice+window,
                     strand = "*"))

fig6region <- as.ggplot(expression(
    kp <- plotKaryotype(genome="susScr11", chromosomes = c("chr2")),
    kpPlotDensity(kp, all.genes, col="grey80", border = "white"),
    kpPlotMarkers(kp, chr=gene_labels$chr, y= 0.5,
                 x=gene_labels$pos, labels=gene_labels$label),
    kpPlotRegions(kp, region)
))

### A: fgsea plot --------------------------------------------------------------

fgsea_plot <- readRDS("data/GG_enrichment_plot.RDS")

fig6c <- ggplot(fgsea_plot$data[[1]]) +
    geom_line(aes(x=x, y=y)) +
    geom_segment(data = fgsea_plot$data[[2]],
                 aes(x=x, xend=xend, y=y, yend=yend), linewidth=0.2, color="grey40") +
    ggrepel::geom_text_repel(data = data.frame(x=142,y=-0.05596333, xend=142,
                                               yend=0.05596333, label=paste0("Variant 2:", pos_choice)),
                             aes(x=x,y=yend, label=label),
                             nudge_x=150000, nudge_y=0.15, color="grey40") +
    geom_hline(aes(yintercept = 0)) +
    theme_minimal() +
    ylab("Enrichment") + xlab("Rank")


### B: Boxplot -----------------------------------------------------------------

BASE <- "WGS_300_GS_filter_chr2"
bim_file <- read.table(paste0("data/", BASE, ".bim"))
fam_file <- read.table(paste0("data/", BASE, ".fam"))
genotypes <- read_bed(bed_file = paste0("data/", BASE, ".bed"),
                      bim_file = paste0("data/", BASE, ".bim"),
                      fam_file = paste0("data/", BASE, ".fam"))
out <- genotypes$x
out[out==2] <- NA
out[out==3] <- 2
X <- matrix(out, ncol=genotypes$n, nrow=genotypes$p, byrow=TRUE)
X <- t(X)
rownames(X) <- fam_file$V2
colnames(X) <- bim_file$V2

gene_choice <- "IGF2"
chr_choice <- "chr2"
tissue_choice <- "liver"
pos <- pos_choice
expression_LW <-
    read.table(paste0(
                      "data/WGS_100.LW_GS_filter_chr2_ENSSSCG00000035293_IGF2_",
                      tissue_choice, ".dat"))
expression_LD <-
    read.table(paste0(
                      "data/WGS_100.LD_GS_filter_chr2_ENSSSCG00000035293_IGF2_",
                      tissue_choice, ".dat"))
expression_DU <-
    read.table(paste0(
                      "data/WGS_100.DU_GS_filter_chr2_ENSSSCG00000035293_IGF2_",
                      tissue_choice, ".dat"))
LW_keep <- read.table(paste0("data/LW.keep"))
LD_keep <- read.table(paste0("data/LD.keep"))
DU_keep <- read.table(paste0("data/DU.keep"))
expression_LW_df <- data.frame(ID = LW_keep$V2, breed=LW_keep$V1,
                               expression=expression_LW$V1)
expression_LD_df <- data.frame(ID = LD_keep$V2, breed=LD_keep$V1,
                               expression=expression_LD$V1)
expression_DU_df <- data.frame(ID = DU_keep$V2, breed=DU_keep$V1,
                               expression=expression_DU$V1)
expression_df <- bind_rows(expression_LW_df, expression_LD_df) %>%
    bind_rows(., expression_DU_df)

index <- which(bim_file$V4 == pos)
geno_choose <- data.frame(genotype = factor(X[,index]), ID = fam_file$V2)
plot_df <- full_join(expression_df, geno_choose, by = "ID")
levels(plot_df$genotype) <- c("AA", "AC", "CC")
fig6d <- ggplot(data = plot_df) +
    geom_boxplot(aes(x=genotype, y=expression, fill=breed), alpha = 0.5) +
    geom_point(aes(x=genotype, y=expression), alpha = 0.2) +
    geom_hline(aes(yintercept = 0), lty = 2) +
    facet_wrap(~breed) +
    theme_bw()+
    guides(fill="none") +
    ylab("Expression") +
    xlab(paste0("Variant 2:", pos, " genotype"))

### C: Track -------------------------------------------------------------------

## Read in data
ensembl <- rtracklayer::import('data/susScr11.ensGene.gtf.gz')
df_ensembl = as.data.frame(ensembl)
df_ensembl$strand = "*"
gr_ensembl = makeGRangesFromDataFrame(df_ensembl, keep.extra.columns = T)
df_cat <- read.table("data/Traits_category.txt",
                     header=T,
                     sep="\t")
QTLdb <- read.table("data/QTLdb_pigSS11.bed",
                    skip = 18,
                    fill = T,
                    row.names = NULL,
                    header = F,
                    sep="\t")
QTLdb = QTLdb[,-c(5:12)]
colnames(QTLdb) = c("chromosome","start","end","trait")
QTLdb$trait = sub("\\(.*", "", QTLdb$trait)
QTLdb$chromosome = sub("....","", QTLdb$chromosome)
QTLdb$full_trait = gsub(" QTL ", "",QTLdb$trait)
QTLdb = left_join(QTLdb,
                  df_cat)
gr_QTLdb = makeGRangesFromDataFrame(QTLdb, keep.extra.columns = T)
gr_QTLdb_filtered = makeGRangesFromDataFrame(filter(as.data.frame(gr_QTLdb), width <10),
                                             keep.extra.columns = T)

ensembl_choice <- ensembl %>% as.data.frame() %>%
    filter(seqnames == chr_choice, type == "transcript") %>%
    group_by(gene_id) %>%
    summarize(min_start = min(start), max_end = max(end)) %>%
    ungroup() %>%
    filter(min_start > pos_choice-window, max_end < pos_choice+window) %>%
    separate(gene_id, into=c("ensembl_gene_id", "misc"))
ensembl_to_genename <- read.table("data/liver_ensembl-to-genename.txt", header=TRUE,
                                  sep="\t", quote="")
ensembl_choice <- left_join(ensembl_choice, ensembl_to_genename,
                            by = "ensembl_gene_id")

QTLdb_choice <- gr_QTLdb_filtered %>% as.data.frame() %>%
    filter(seqnames == 2, start > pos_choice-window,
           end < pos_choice+window) %>%
    dplyr::select(start, end, cat_trait) %>% unique()
QTLdb_choice$cat_trait <- factor(QTLdb_choice$cat_trait)
levels(QTLdb_choice$cat_trait) <- c("Health QTL", "Meat and carcass QTL")

## chr2:50750725  (DU) => fort Vbeta and QTL annotated
start <- pos_choice
end <- pos_choice

full_annot <-
    read.table(paste0("data/",
                      "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
                      tissue_choice, "-stages123_", chr_choice, ".txt"))
colnames(full_annot) <-
    c("REG", "NSC", "OCR_30dpf", "OCR_70dpf", "OCR_newborn",
      "UMR_30dpf", "LMR_30dpf", "UMR_70dpf", "LMR_70dpf",
      "UMR_newborn", "LMR_newborn", "other")
annot_names <-
    read.table(paste0("data/",
                      "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
                      tissue_choice, "-stages123_", chr_choice,
                      "_annotation-order.txt"))

pos <- read.table(paste0("data/WGS_300_GS_filter_", chr_choice, ".bim"))

Vbeta_DU <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_DU.param"), header=TRUE)
Vbeta_LD <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_LD.param"), header=TRUE)
Vbeta_LW <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_LW.param"), header=TRUE)
Vbeta_DU <- data.frame(pos=pos$V4, Vbeta=Vbeta_DU$Vbeta, breed="DU")
Vbeta_LD <- data.frame(pos=pos$V4, Vbeta=Vbeta_LD$Vbeta, breed="LD")
Vbeta_LW <- data.frame(pos=pos$V4, Vbeta=Vbeta_LW$Vbeta, breed="LW")
Vbeta <- bind_rows(Vbeta_DU, Vbeta_LD) %>% bind_rows(., Vbeta_LW)

## Make plot
full_annot_long <- gather(cbind(full_annot[,-which(colnames(full_annot) == "other")],
                                pos=pos$V4),
                          key = "stage", value = "annotation", -pos)
datfilter <- full_annot_long %>%
    filter(annotation != 0) %>%
    filter(pos > (start-window), pos < (end+window))
datfilter$stage <- factor(datfilter$stage)
datfilter$y <- as.numeric(datfilter$stage)
datfilter <- rbind(datfilter,data.frame(pos = c(start,end),
                                        stage=gene_choice, annotation=1, y=0))
datfilter$stage <- fct_relevel(datfilter$stage,
                               "REG", "NSC",
                               "OCR_30dpf", "OCR_70dpf", "OCR_newborn",
                               "UMR_30dpf", "UMR_70dpf", "UMR_newborn",
                               "LMR_30dpf", "LMR_70dpf", "LMR_newborn")
levels(datfilter$stage) <- gsub("_", " ", levels(datfilter$stage))
fig6a0 <- ggplot(datfilter %>% filter(stage != gene_choice)) +
    annotate("rect", fill="grey90",  xmin=start, xmax=end,
             ymin = -Inf,ymax=Inf) +
    geom_point(data=QTLdb_choice %>% mutate(mean_pos=(start+end)/2),
               aes(x=mean_pos, y=factor(cat_trait))) +
    geom_vline(aes(xintercept=pos_choice), linetype=2) +
    ylab("") +
    scale_color_viridis(discrete=TRUE,direction = -1) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    guides(color = "none") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
fig6a1 <- ggplot(datfilter %>% filter(stage != gene_choice)) +
    annotate("rect", fill="grey90",  xmin=start, xmax=end,
             ymin = -Inf,ymax=Inf) +
    geom_point(aes(y=stage, x=pos)) +
    geom_vline(aes(xintercept=pos_choice), linetype=2) +
    ylab("") +
    xlab(paste0("Position on chromosome ", substr(chr_choice, 4,5))) +
    scale_color_viridis(discrete=TRUE,direction = -1) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    guides(color = "none") +
    theme_bw()
fig6a2 <- ggplot(Vbeta %>% filter(pos > (start-window), pos < (end+window))) +
    annotate("rect", fill="grey90",  xmin=start, xmax=end,
             ymin = -Inf,ymax=Inf) +
    geom_point(data = Vbeta %>%
                   filter(pos > (start-window), pos < (end+window)),
               aes(y=Vbeta, x=pos, color=breed), alpha = 0.5) +
    # geom_vline(data = Vbeta %>%
    #                filter(pos > (start-window), pos < (end+window), eQTL == 1),
    #            aes(xintercept=pos), alpha = 0.5) +
    facet_wrap(~breed, ncol=1, scales="free_y", strip.position = "right") +
    ylab("Posterior variance") + xlab("") +
    # scale_x_discrete(labels=levels(datfilter$stage)) +
    guides(color="none") +
    geom_vline(aes(xintercept=pos_choice), linetype=2) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    theme_bw() +
    theme(strip.background = element_blank(), axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

fig6a3 <- ggplot(datfilter %>% filter(stage != gene_choice)) +
    geom_segment(data = ensembl_choice,
                 aes(x=min_start, xend=max_end, y=0, yend=0), linewidth = 4) +
    ggrepel::geom_text_repel(data = ensembl_choice %>%
                                 mutate(middle = (min_start+max_end)/2),
                             aes(x=middle, y=0, label=hgnc_symbol), nudge_y=0.5,
                             force=2, size = 3) +
    ylab("") +
    scale_color_viridis(discrete=TRUE,direction = -1) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    guides(color = "none") +
    theme_nothing() +
    theme(strip.background = element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
## 816 x 364
fig6a <- plot_grid(#fig6a3,
                   fig6a2, fig6a0, fig6a1,
                   ncol=1, align = "v", axis = 'lr',
                   rel_heights = c(2, 0.5, 2))

### D: Medium/high class -------------------------------------------------------

annot_file <- read.table("data/allchr_annotations.txt",header=T)
annot_file$Annotation = rep("Not annotated", nrow(annot_file))
annot_file$Annotation[which(annot_file$REG == 1 | annot_file$NSC ==1)] = "REG+NSC"
annot_file$Annotation[which(annot_file$REG == 0 & annot_file$NSC == 0 &
                                annot_file$other == 0 &
                                annot_file$tissue == "liver")]= "liver-specific"
annot_file$Annotation[which(annot_file$REG == 0 & annot_file$NSC == 0 &
                                annot_file$other == 0 &
                                annot_file$tissue == "muscle")]= "muscle-specific"

IGF2_liver_LD_param = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_LD.param",
                                 header=T)
IGF2_liver_LW_param = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_LW.param",
                                 header=T)
IGF2_liver_DU_param = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_DU.param",
                                 header=T)

IGF2_liver_LD_hyp = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_LD.hyp",
                               header=T)
IGF2_liver_LW_hyp = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_LW.hyp",
                               header=T)
IGF2_liver_DU_hyp = 
  read.table("data/chr2_IGF2_liver_BayesRCO_vep.atacseq.methylation3.matched_DU.hyp",
                               header=T)

df_IGF2_liver_LD = hyp_to_annot_enrichment(IGF2_liver_LD_hyp, nclass= 5, 
                                           colnames(annot_file)[3:14])
df_IGF2_liver_LD$Class = as.factor(df_IGF2_liver_LD$Class)
df_IGF2_liver_LD$Annotation = 
  factor(df_IGF2_liver_LD$Annotation,
         levels=arrange(aggregate(filter(df_IGF2_liver_LD, 
                                         Class %in% c(3,4,5))$Freq, 
                                  list(filter(df_IGF2_liver_LD, 
                                              Class %in% c(3,4,5))$Annotation),sum),x)$Group.1)
df_IGF2_liver_LW = 
  hyp_to_annot_enrichment(IGF2_liver_LW_hyp, nclass= 5, colnames(annot_file)[3:14])
df_IGF2_liver_DU = 
  hyp_to_annot_enrichment(IGF2_liver_DU_hyp, nclass= 5, colnames(annot_file)[3:14])

df_IGF_liver = rbind(cbind(df_IGF2_liver_DU, Breed="DU"),
                     cbind(df_IGF2_liver_LW, Breed="LW"),
                     cbind(df_IGF2_liver_LD, Breed="LD"))
df_IGF = cbind(df_IGF_liver, Tissue ="liver")

df_annot_rename = data.frame(Annotation=unique(df_IGF$Annotation),
                             Annotation2= c("REG","NSC",
                                            "Peaks 30dpf","Peaks 70dpf", "Peaks NB",
                                            "UMR 30dpf", "LMR 30dpf",
                                            "UMR 70dpf", "LMR 70dpf",
                                            "UMR NB", "LMR NB",
                                            "Unannotated"))
df_IGF = left_join(df_IGF, df_annot_rename)
df_IGF$Annotation2 = 
  factor(df_IGF$Annotation2,
         levels=df_annot_rename[match(arrange(aggregate(filter(df_IGF, Breed == "DU" & Tissue =="liver" &  Class %in% c(3,4,5))$Freq, 
                                                        list(filter(df_IGF, Breed == "DU" & Tissue =="liver" &  Class %in% c(3,4,5))$Annotation),sum),x)$Group.1,
                                      df_annot_rename$Annotation),]$Annotation2)

df_class = data.frame(Class = as.character(c(1,2,3,4,5)),
                      ClassSize = c("Null","Very low", "Low", "Medium", "High"))

df_IGF = left_join(df_IGF,df_class)
df_IGF$ClassSize = factor(df_IGF$ClassSize,
                          levels=c("Low","Medium","High"))
df_IGF$Annotation2 = factor(df_IGF$Annotation2,
                            levels=df_annot_rename[match(arrange(aggregate(filter(df_IGF, Breed == "DU" & Tissue =="liver" & Class %in% c(4,5))$Freq, list(filter(df_IGF, Breed == "DU" & Tissue =="liver" &  Class %in% c(4,5))$Annotation),sum),x)$Group.1,df_annot_rename$Annotation),]$Annotation2)
levels(df_IGF$Annotation2)[c(6,7,9,10,11)] <-
    c("OCR newborn", "LMR newborn", "OCR 30dpf", "OCR 70dpf", "UMR newborn")
fig6b <- ggplot(filter(df_IGF, Class %in% c(4,5), Tissue == "liver"),
       aes(x=Annotation2, y=Freq, fill=Breed, alpha=ClassSize))+
    geom_bar(stat="identity")+
    theme_bw()+
    coord_flip()+
    # facet_rep_grid(Tissue~Breed)+
    facet_rep_grid(~Breed)+
    scale_alpha_discrete(range = c(0.5,1)) +
    theme(strip.background = element_blank(),
          #      panel.grid.minor = element_blank(),
          #     panel.grid.major = element_blank(),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=12),
          strip.text = element_text(size=12))+
    xlab("")+
    ylab("Frequency of high and medium effect variants")+
    theme(legend.position="none") +
    guides(scale="none")
    # scale_fill_manual(name = "Effect size", values=c("grey70", "black"))

### Put it all together --------------------------------------------------------
fig6bot<- plot_grid(fig6a, fig6b, labels = c("C","D"), ncol=1,
                     rel_heights=c(2.2,1))
fig6top <- plot_grid(fig6c, fig6d, labels = c("A", "B"), nrow = 1)
fig6ttop <- plot_grid(fig6region, fig6top, nrow = 2,
                      rel_heights = c(.5,1))
## 750 x 1100
fig6 <- plot_grid(fig6ttop, fig6bot, nrow = 2, rel_heights = c(1,3))

## Supplementary Figure --------------------------------------------------------

### Supp Figure 1 (Density of SNPs after filtering  ----------------------------
bim_file <- read.table(paste0("data/WGS_300_GS_filter_chr", 1, ".bim"))
for(i in c(2,5,6,7,10,14,18)) {
  tmp <- read.table(paste0("data/WGS_300_GS_filter_chr", i, ".bim"))
  bim_file <- bind_rows(bim_file, tmp)
}
bim_file$chr <- factor(paste0("chr", bim_file$V1))
bim_file$chr <- fct_relevel(bim_file$chr,
                            paste0("chr", c(1,2,5,6,7,10,14,18)))
ggplot(bim_file) +
  geom_density(aes(x=V4)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  xlab("Position") +
  ylab("Density") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Supp Figure 2 (Per-chr genomic PCAs) ---------------------------------------
eigenvec <- vector("list", length=8)
names(eigenvec) <- paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18))
eigenval <- vector("list", length=8)
names(eigenval) <- paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18))
for(chr in names(eigenvec)) {
  eigenvec[[chr]] <- read.table(paste0("data/WGS_300_GS_filter_prune_",
                                       chr, ".eigenvec"))
  eigenval[[chr]] <- read.table(paste0("data/WGS_300_GS_filter_prune_",
                                       chr, ".eigenval"))
}

eigenvec_df <- bind_rows(eigenvec, .id="chr")
colnames(eigenvec_df)[c(2:5)] <- c("breed", "ID", "PCA1", "PCA2")
eigenvec_df$chr <- factor(eigenvec_df$chr,
                          levels = paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18)))

perc1 <- lapply(eigenval, function(x) {
  data.frame(PC1 = paste0("PC1: ",round(x$V1[1] / sum(x$V1)*100,1), "%"))}) %>%
  bind_rows(., .id="chr")
perc2 <- lapply(eigenval, function(x) {
  data.frame(PC2 = paste0("PC2: ",round(x$V1[2] / sum(x$V1)*100,1), "%"))}) %>%
  bind_rows(., .id="chr")
perc <- full_join(perc1, perc2, by="chr")
perc$perc <- paste0(perc$PC1, "\n", perc$PC2)
perc$chr <- factor(perc$chr,
                   levels = paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18)))

ggplot(eigenvec_df) +
  geom_point(aes(x=PCA1, y=PCA2, col=breed), alpha = 0.3) +
  geom_hline(aes(yintercept = 0), lty=2, col="grey30") +
  geom_vline(aes(xintercept = 0), lty=2, col="grey30") +
  geom_text(data=perc, aes(x=Inf, y=-Inf, label=perc), hjust = 1, vjust = -.1) +
  facet_wrap(~chr) +
  theme_bw() +
  scale_color_discrete(name="") +
  theme(legend.position="bottom")


### Supp Figure 3 (sex-corrected expression) -----------------------------------
all.files <-
  data.frame(file = list.files("data/expression_data/")) %>%
  mutate(gene = split_choose(file, split="_", pos = 7),
         tissue = split_choose(file, split="_", pos = 8),
         pop = substr(split_choose(file, split="_", pos = 2), 5,20)) %>%
  separate(tissue, into=c("tissue", "misc")) %>%
  dplyr::select(-misc) %>%
  filter(pop %in% c("DU", "LD", "LW"))

all.dat <- data.frame(gene = character(), tissue = character(),
                      pop = character(), expression = numeric())
for(g in unique(all.files$gene)) {
  for(ti in unique(all.files$tissue)) {
    tmp <- all.files %>% filter(gene == g, tissue == ti)
    if(nrow(tmp) == 0) next;
    for(i in 1:nrow(tmp)) {
      tmp2 <- read.table(paste0("data/expression_data/", tmp$file[i]))$V1
      all.dat <- rbind(all.dat, data.frame(gene = g, tissue = ti,
                                           pop = tmp$pop[i],
                                           expression = tmp2))
    }
  }
}

all.dat$pop <- factor(all.dat$pop)
levels(all.dat$pop) <- c("Duroc", "Landrace", "LargeWhite")

## Remove L3HYPDH because no longer included in Crespo et al
g <- ggplot(all.dat %>% filter(gene != "L3HYPDH")) +
  geom_boxplot(aes(x=gene, y=expression, fill = pop, color= pop), alpha = 0.3) +
  geom_hline(aes(yintercept=0), lty = 2, col="black") +
  facet_wrap(~tissue, ncol=1) +
  theme_bw() +
  ylab("Sex-corrected logCPM") + xlab("") +
  theme(legend.position="none")


### Supp Figure 4 (per-chr barplots of annotations) -----------------------------

dat <- bind_rows(data.frame(df_left_muscle, tissue = "muscle"),
                 data.frame(df_left_liver, tissue = "liver")) %>%
  filter(QTL == "Yes")

dat <- dat %>%
  dplyr::select(Chromosome = seqnames, `Trait category`=cat_trait.y, Tissue = tissue.x,
         Annotation = Annotation.x, QTL = QTL)

dat$`Trait category` <- factor(dat$`Trait category`)
levels(dat$`Trait category`)[3] <- "Meat and Carcass"

plotdat <- dat %>%
  group_by(Chromosome, Tissue, `Trait category`) %>%
  summarize(n=n()) %>%
  ungroup()
plotdat$Chromosome <- fct_rev(factor(plotdat$Chromosome))

ggplot(plotdat) +
  geom_col(aes(x=Chromosome, y=n, fill=`Trait category`),
           position = position_dodge()) +
  ylab("Number of QTLs") +
  facet_grid(~Tissue) +
  coord_flip() +
  theme_bw() +
  theme(legend.position="bottom") 


### Supp Figure 5 (Fisher's exact test) ----------------------------------------

annot_pos = data.frame(seqnames = annot_bim$chr,
                       start = annot_bim$bp_pos,
                       end = annot_bim$bp_pos,
                       strand = NA,
                       allele_1 = annot_bim$allele_1,
                       allele_2 =  annot_bim$allele_2)

df_annot_muscle = cbind(annot_pos,
                        filter(annot_file, tissue=="muscle"))
df_annot_liver = cbind(annot_pos,
                       filter(annot_file, tissue=="liver"))

gr_annot_muscle = makeGRangesFromDataFrame(df_annot_muscle,
                                           keep.extra.columns = T,
                                           seqnames.field="seqnames")
gr_annot_liver = makeGRangesFromDataFrame(df_annot_liver,
                                          keep.extra.columns = T,
                                          seqnames.field="seqnames")

gr_full_muscle = join_overlap_full(gr_annot_muscle,gr_QTLdb_filtered)
gr_full_filtered_muscle = subsetByOverlaps(subsetByOverlaps(gr_full_muscle, 
                                                            gr_annot_muscle), gr_QTLdb_filtered)

gr_left_muscle = subsetByOverlaps(gr_full_muscle, gr_annot_muscle)
df_left_muscle =as.data.frame(gr_left_muscle)
df_left_muscle$cat_trait.y[is.na(df_left_muscle$cat_trait.y)]="None"
contigency_table = table(df_left_muscle$cat_trait.y, df_left_muscle$Annotation.x)

df_fisher_muscle = data.frame(Category = rownames(contigency_table),
                              pval = rep(1, 6),
                              OR = rep(0,6),
                              IC_min = rep(0,6),
                              IC_max = rep(0,6))
df_fisher_muscle$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$p.value})
df_fisher_muscle$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$estimate})
df_fisher_muscle$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[1]})
df_fisher_muscle$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[2]})


df_QTL_muscle = as.data.frame(gr_full_filtered_muscle)
df_QTL_muscle$Freq = rep(1, nrow(df_QTL_muscle))
freq_QTL_annot_muscle = aggregate(df_QTL_muscle$Freq,
                                  list(df_QTL_muscle$Annotation.x, df_QTL_muscle$cat_trait.y),
                                  sum)
colnames(freq_QTL_annot_muscle) = c("Annotated","Trait category","Total")
freq_QTL_annot_muscle$Freq = freq_QTL_annot_muscle$Total / rep(table(df_annot_muscle$Annotation),5)


gr_full_liver = join_overlap_full(gr_annot_liver,gr_QTLdb_filtered)
gr_full_filtered_liver = subsetByOverlaps(subsetByOverlaps(gr_full_liver, gr_annot_liver), gr_QTLdb_filtered)

gr_left_liver = subsetByOverlaps(gr_full_liver, gr_annot_liver)
df_left_liver =as.data.frame(gr_left_liver)
df_left_liver$cat_trait.y[is.na(df_left_liver$cat_trait.y)]="None"
contigency_table = table(df_left_liver$cat_trait.y, df_left_liver$Annotation.x)

df_fisher_liver = data.frame(Category = rownames(contigency_table),
                             pval = rep(1, 6),
                             OR = rep(0,6),
                             IC_min = rep(0,6),
                             IC_max = rep(0,6))
df_fisher_liver$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$p.value})
df_fisher_liver$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$estimate})
df_fisher_liver$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[1]})
df_fisher_liver$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[2]})

df_fisher_REG = data.frame(Category = rownames(contigency_table),
                           pval = rep(1, 6),
                           OR = rep(0,6),
                           IC_min = rep(0,6),
                           IC_max = rep(0,6))
df_fisher_REG$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$p.value})
df_fisher_REG$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$estimate})
df_fisher_REG$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$conf.int[1]})
df_fisher_REG$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$conf.int[2]})

df_fisher_all = rbind(cbind(df_fisher_muscle, Tissue = "Muscle"),
                      cbind(df_fisher_liver, Tissue = "Liver"),
                      cbind(df_fisher_REG, Tissue = "Not context-specific"))

df_fisher_all$Threshold = "False"
df_fisher_all$Threshold[which(df_fisher_all$pval<0.05)]="True"

df_fisher_all$Category <- factor(df_fisher_all$Category)
levels(df_fisher_all$Category)[3] <- "Meat and Carcass"
ggplot(filter(df_fisher_all,Category != "None"), aes(y=Category, x=OR, label=Category, color = -log10(pval))) +
  geom_point(size=4) +
  geom_errorbarh(aes(xmin=IC_min, xmax=IC_max), height=.3) +
  coord_fixed(ratio=.3) +
  ylab("") +
  geom_vline(xintercept=1, linetype='longdash') +
  scale_color_gradient(low = "#FCBBA1", high="#A50F15", name = "-log10 p-value")+
  facet_wrap(~Tissue, ncol=3)+
  theme_bw()


### Supp Figure 6 (correlation of per-variant V_i) -----------------------------

Vbeta_annot <- read.table("data/Vbeta_annot_all.txt",
                          header=T,
                          sep="\t")
Vbeta_annot = Vbeta_annot[-grep("tissue",Vbeta_annot$tissue),]
Vbeta_annot$start = as.numeric(Vbeta_annot$start)
Vbeta_annot$end = as.numeric(Vbeta_annot$end)
Vbeta_annot$Sum_vi = as.numeric(Vbeta_annot$Sum_vi)
Vbeta_annot$Vbeta = as.numeric(Vbeta_annot$Vbeta)

df_cor_Vbeta = data.frame(Scenario = rep(NA, 66),
                          Cor_Vbeta = rep(NA, 66),
                          Cor_Vbeta_sum = rep(NA, 66))
compt = 1
for (i in unique(Vbeta_annot$gene)){
  for (j in unique(Vbeta_annot$tissue)){
    temp = filter(Vbeta_annot, gene== i & tissue==j)
    df_cor_Vbeta$Scenario[compt] = paste0(i, "_", j, "_DU_LD",sep="")
    df_cor_Vbeta$Scenario[compt+1] = paste0(i, "_", j, "_DU_LW",sep="")
    df_cor_Vbeta$Scenario[compt+2] = paste0(i, "_", j, "_LW_LD",sep="")
    df_cor_Vbeta$Cor_Vbeta_sum[compt] =
      cor(as.numeric(filter(temp, learning=="DU")$Sum_vi),
          as.numeric(filter(temp, learning=="LD")$Sum_vi))
    df_cor_Vbeta$Cor_Vbeta_sum[compt+1] = 
      cor(as.numeric(filter(temp, learning=="DU")$Sum_vi),
          as.numeric(filter(temp, learning=="LW")$Sum_vi))
    df_cor_Vbeta$Cor_Vbeta_sum[compt+2] = 
      cor(as.numeric(filter(temp, learning=="LW")$Sum_vi),
          as.numeric(filter(temp, learning=="LD")$Sum_vi))
    df_cor_Vbeta$Cor_Vbeta[compt] = 
      cor(as.numeric(filter(temp, learning=="DU")$Vbeta),
          as.numeric(filter(temp, learning=="LD")$Vbeta))
    df_cor_Vbeta$Cor_Vbeta[compt+1] = 
      cor(as.numeric(filter(temp, learning=="DU")$Vbeta),
          as.numeric(filter(temp, learning=="LW")$Vbeta))
    df_cor_Vbeta$Cor_Vbeta[compt+2] = 
      cor(as.numeric(filter(temp, learning=="LW")$Vbeta),
          as.numeric(filter(temp, learning=="LD")$Vbeta))
    compt = compt + 3
  }
}

df_cor_Vbeta$Gene = rep(unique(Vbeta_annot$gene), each = 6)
df_cor_Vbeta$Tissue = rep(unique(Vbeta_annot$tissue), each = 3)
df_cor_Vbeta$Breed1 = rep(c("DU","DU","LW"))
df_cor_Vbeta$Breed2 = rep(c("LD","LW","LD"))
df_cor_Vbeta$Pairwise = rep(c("DU_LD","DU_LW","LW_LD"))
df_cor_Vbeta$Scenario1 = paste0(df_cor_Vbeta$Gene,
                                "_",
                                df_cor_Vbeta$Tissue,
                                "_",
                                df_cor_Vbeta$Breed1,
                                sep="")
df_cor_Vbeta$Scenario2 = paste0(df_cor_Vbeta$Gene,
                                "_",
                                df_cor_Vbeta$Tissue,
                                "_",
                                df_cor_Vbeta$Breed2,
                                sep="")
df_cor_Vbeta$Pairwise <- factor(df_cor_Vbeta$Pairwise)
levels(df_cor_Vbeta$Pairwise) <- c("(DU,LD)", "(DU,LW)", "(LW,LD)")

## Remove L3HYPDH because no longer included in Crespo et al
supfig6 <- ggplot(df_cor_Vbeta %>% filter(Gene != "L3HYPDH"),
                aes(x=Cor_Vbeta, color = Pairwise, fill=Pairwise))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Tissue)+
  xlab(expression(paste("Correlation between individual ", V[i])))+
  xlim(c(0,1))+
  scale_colour_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"),
                      name="")+
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"),
                    name="") +
  theme(legend.position="bottom")




### Supp Figure 7+8 (Summed window-based V_values) -----------------------------

window_vi <- read.table("data/Vbeta_annot_sum_IGF2.txt",
                        header=TRUE)
setDT(window_vi)[,groupid := rleid(Sum_vi)][]
window_vi_df <- as.data.frame(window_vi)

window_vi_agg <- window_vi_df %>%
  group_by(tissue, learning, groupid) %>%
  summarise(Sum_Vi_unique = unique(Sum_vi), start = mean(start))


gene_choice <- "IGF2"
ensembl_choice <- "ENSSSCG00000035293"
chr_choice <- "chr2"
tissue_choice <- "liver"
ens <- read.table("data/susScr11.ensGene.gtf", header=FALSE,
                  sep="\t")
exons <- ens %>% filter(grepl(ensembl_choice, V9), V3 == "exon") %>%
  dplyr::select(V4, V5) %>% unique()
gene <- ens %>% filter(grepl(ensembl_choice, V9)) %>%
  dplyr::select(V4, V5) %>% unique()
start <- min(gene$V4)
end <- max(gene$V5)
middle <- mean(c(start, end))


figsxa <- ggplot(data = window_vi_agg %>% filter(tissue == "liver")) +
  geom_point(aes(x=start, y=Sum_Vi_unique, col=learning), alpha = 0.3) +
  geom_vline(aes(xintercept = middle), linetype = 2) +
  facet_wrap(~learning, ncol=1, scales = "free_y") +
  xlab("Position on chromosome 2") +
  ylab("Window-based posterior variance estimate") +
  guides(color="none") +
  theme_bw()


figsxb <- ggplot(data = window_vi_agg %>% filter(tissue == "muscle")) +
  geom_point(aes(x=start, y=Sum_Vi_unique, col=learning), alpha = 0.3) +
  geom_vline(aes(xintercept = middle), linetype = 2) +
  facet_wrap(~learning, ncol=1, scales = "free_y") +
  xlab("Position on chromosome 2") +
  ylab("Window-based posterior variance estimate") +
  scale_color_discrete(name = "") +
  guides(color="none") +
  theme_bw()


### Supp Figure 9 (IGF2 track) -------------------------------------------------
gene_choice <- "IGF2"
ensembl_choice <- "ENSSSCG00000035293"
chr_choice <- "chr2"
tissue_choice <- "liver"
window <- 1000000 #100000

## Read in data
ensembl <- rtracklayer::import('data/susScr11.ensGene.gtf.gz')
df_ensembl = as.data.frame(ensembl)
df_ensembl$strand = "*"
gr_ensembl = makeGRangesFromDataFrame(df_ensembl, keep.extra.columns = T)
df_cat <- read.table("data/Traits_category.txt",
                     header=T,
                     sep="\t")
QTLdb <- read.table("data/QTLdb_pigSS11.bed",
                    skip = 18,
                    fill = T,
                    row.names = NULL,
                    header = F,
                    sep="\t")
QTLdb = QTLdb[,-c(5:12)]
colnames(QTLdb) = c("chromosome","start","end","trait")
QTLdb$trait = sub("\\(.*", "", QTLdb$trait)
QTLdb$chromosome = sub("....","", QTLdb$chromosome)
QTLdb$full_trait = gsub(" QTL ", "",QTLdb$trait)
QTLdb = left_join(QTLdb,
                  df_cat)
gr_QTLdb = makeGRangesFromDataFrame(QTLdb, keep.extra.columns = T)
gr_QTLdb_filtered = makeGRangesFromDataFrame(filter(as.data.frame(gr_QTLdb), width <10),
                                             keep.extra.columns = T)

pos_choice <- mean(c(1469132,1496442)) 


ensembl_choice <- ensembl %>% as.data.frame() %>%
    filter(seqnames == chr_choice, type == "transcript") %>%
    group_by(gene_id) %>%
    summarize(min_start = min(start), max_end = max(end)) %>%
    ungroup() %>%
    filter(min_start > pos_choice-window, max_end < pos_choice+window) %>%
    separate(gene_id, into=c("ensembl_gene_id", "misc"))
ensembl_to_genename <- read.table("data/liver_ensembl-to-genename.txt", header=TRUE,
                                  sep="\t", quote="")
ensembl_choice <- left_join(ensembl_choice, ensembl_to_genename,
                            by = "ensembl_gene_id")
ensembl_choice$hgnc_symbol <- ifelse(ensembl_choice$wikigene_name == "IGF2",
                                     "IGF2", ensembl_choice$hgnc_symbol)

QTLdb_choice <- gr_QTLdb_filtered %>% as.data.frame() %>%
    filter(seqnames == 2, start > pos_choice-window,
           end < pos_choice+window) %>%
    dplyr::select(start, end, cat_trait) %>% unique()
QTLdb_choice$cat_trait <- factor(QTLdb_choice$cat_trait)
levels(QTLdb_choice$cat_trait) <- c("Meat and carcass QTL")

start <- pos_choice
end <- pos_choice

full_annot <-
    read.table(paste0("data/",
                      "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
                      tissue_choice, "-stages123_", chr_choice, ".txt"))
colnames(full_annot) <-
    c("REG", "NSC", "OCR_30dpf", "OCR_70dpf", "OCR_newborn",
      "UMR_30dpf", "LMR_30dpf", "UMR_70dpf", "LMR_70dpf",
      "UMR_newborn", "LMR_newborn", "other")
annot_names <-
    read.table(paste0("data/",
                      "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
                      tissue_choice, "-stages123_", chr_choice,
                      "_annotation-order.txt"))

pos <- read.table(paste0("data/WGS_300_GS_filter_", chr_choice, ".bim"))

Vbeta_DU <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_DU.param"), header=TRUE)
Vbeta_LD <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_LD.param"), header=TRUE)
Vbeta_LW <-
    read.table(paste0("data/", chr_choice, "_", gene_choice, "_",
                      tissue_choice, "_BayesRCO_vep.atacseq.methylation3.matched",
                      "_LW.param"), header=TRUE)
Vbeta_DU <- data.frame(pos=pos$V4, Vbeta=Vbeta_DU$Vbeta, breed="DU")
Vbeta_LD <- data.frame(pos=pos$V4, Vbeta=Vbeta_LD$Vbeta, breed="LD")
Vbeta_LW <- data.frame(pos=pos$V4, Vbeta=Vbeta_LW$Vbeta, breed="LW")
Vbeta <- bind_rows(Vbeta_DU, Vbeta_LD) %>% bind_rows(., Vbeta_LW)

## Make plot
full_annot_long <- gather(cbind(full_annot[,-which(colnames(full_annot) == "other")],
                                pos=pos$V4),
                          key = "stage", value = "annotation", -pos)
datfilter <- full_annot_long %>%
    filter(annotation != 0) %>%
    filter(pos > (start-window), pos < (end+window))
datfilter$stage <- factor(datfilter$stage)
datfilter$y <- as.numeric(datfilter$stage)
datfilter <- rbind(datfilter,data.frame(pos = c(start,end),
                                        stage=gene_choice, annotation=1, y=0))
datfilter$stage <- fct_relevel(datfilter$stage,
                               "REG", "NSC",
                               "OCR_30dpf", "OCR_70dpf", "OCR_newborn",
                               "UMR_30dpf", "UMR_70dpf", "UMR_newborn",
                               "LMR_30dpf", "LMR_70dpf", "LMR_newborn")
levels(datfilter$stage) <- gsub("_", " ", levels(datfilter$stage))
fig6ssa1 <- ggplot(datfilter %>% filter(stage != gene_choice)) +
    annotate("rect", fill="grey90",  xmin=start, xmax=end,
             ymin = -Inf,ymax=Inf) +
    # geom_point(aes(y=stage, x=pos, color = stage)) +
    geom_point(aes(y=stage, x=pos)) +
    geom_point(data=QTLdb_choice %>% mutate(mean_pos = (start+end)/2),
               aes(x=mean_pos, y=factor(cat_trait))) +
    # geom_segment(data = ensembl_choice,
    #              aes(x=min_start, xend=max_end, y=0, yend=0), linewidth = 4) +
    # ggrepel::geom_text_repel(data = ensembl_choice %>%
    #                               mutate(middle = (min_start+max_end)/2),
    #                           aes(x=middle, y=0, label=hgnc_symbol),
    #                          max.overlaps=100, nudge_y=1) +
    ylab("") +
    xlab(paste0("Position on chromosome ", substr(chr_choice, 4,5))) +
    scale_color_viridis(discrete=TRUE,direction = -1) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    guides(color = "none") +
    theme_bw()

fig6ssa2 <- ggplot(Vbeta %>% filter(pos > (start-window), pos < (end+window))) +
    annotate("rect", fill="grey90",  xmin=start, xmax=end,
             ymin = -Inf,ymax=Inf) +
    geom_point(data = Vbeta %>%
                   filter(pos > (start-window), pos < (end+window)),
               aes(y=Vbeta, x=pos, color=breed), alpha = 0.5) +
    # geom_vline(data = Vbeta %>%
    #                filter(pos > (start-window), pos < (end+window), eQTL == 1),
    #            aes(xintercept=pos), alpha = 0.5) +
    facet_wrap(~breed, ncol=1, scales="free_y", strip.position = "right") +
    ylab("Posterior variance") + xlab("") +
    # scale_x_discrete(labels=levels(datfilter$stage)) +
    guides(color="none") +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    theme_bw() +
    theme(strip.background = element_blank(), axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

fig6ssa3 <- ggplot(datfilter %>% filter(stage != gene_choice)) +
    geom_segment(data = ensembl_choice,
                 aes(x=min_start, xend=max_end, y=0, yend=0), linewidth = 4) +
    ggrepel::geom_text_repel(data = ensembl_choice %>%
                                 mutate(middle = (min_start+max_end)/2),
                             aes(x=middle, y=0, label=hgnc_symbol), nudge_y=0.5,
                             force=2, size = 3) +
    ylab("") +
    scale_color_viridis(discrete=TRUE,direction = -1) +
    xlim(c(min(datfilter$pos), max(datfilter$pos))) +
    guides(color = "none") +
    theme_nothing() +
    theme(strip.background = element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
## 816 x 364
fig6ssa <- plot_grid(fig6ssa3, fig6ssa2, fig6ssa1, ncol=1, align = "v", axis = 'lr',
                   rel_heights = c(0.5, 2, 2))








## Supplementary Tables --------------------------------------------------------

### Number of annotated variants -----------------------------------------------

# REG
# NSC
# liver_stage1
# liver_stage2
# liver_stage3
# liver_30dpf_UMR
# liver_30dpf_LMR
# liver_70dpf_UMR
# liver_70dpf_LMR
# liver_NB_UMR
# liver_NB_LMR
# other
res <- data.frame(chr=numeric(), number_variants=numeric(),
                  ocr=character(),
                  umr=numeric(),
                  lmr=numeric(),
                  vep=numeric(),
                  tissue=character())
index <- 1
for(chr in c(1, 2, 5, 6, 7, 10, 14, 18)) {
    for(tissue in c("liver", "muscle")) {
        tmp <- read.table(paste0("WP2_annotations/vep.atacseq.methylation3.matched/",
                                 "vep.atacseq.methylation3.matched_annot_WGS_GS_filter_",
                                 tissue, "-stages123_chr", chr, ".txt"))
        res[index,] <- c(chr, nrow(tmp),
                         length(which(rowSums(tmp[,3:5])>0)),
                         length(which(rowSums(tmp[,c(6,8,10)])>0)),
                         length(which(rowSums(tmp[,c(7,9,11)])>0)),
                         length(which(rowSums(tmp[,1:2])>0)),
                         tissue)

        index <- index + 1
    }
}

### GSEA enrichments -----------------------------------------------------------

res <- read.csv("C:/Users/araul/Desktop/GS_data/fgsea_results/Data/df_bayes_models_list_fgsea_scenarios_padj0.1.csv")

res_1 <- res %>% filter(model == "BayesRCO") %>%
    dplyr::select(-scenario, -pval, -model) %>%
    mutate(padj = round(padj, 3),
           log2err = round(log2err, 3),
           ES = round(ES,3),
           NES = round(NES, 3)) %>%
    arrange(padj)

res_2 <- res %>% filter(model == "BayesR") %>%
    dplyr::select(-scenario, -pval, -model) %>%
    mutate(padj = round(padj, 3),
           log2err = round(log2err, 3),
           ES = round(ES,3),
           NES = round(NES, 3)) %>%
    arrange(padj)

### QTL Category enrichments ---------------------------------------------------

annot_pos = data.frame(seqnames = annot_bim$chr,
                       start = annot_bim$bp_pos,
                       end = annot_bim$bp_pos,
                       strand = NA,
                       allele_1 = annot_bim$allele_1,
                       allele_2 =  annot_bim$allele_2)

df_annot_muscle = cbind(annot_pos,
                        filter(annot_file, tissue=="muscle"))
df_annot_liver = cbind(annot_pos,
                       filter(annot_file, tissue=="liver"))

gr_annot_muscle = makeGRangesFromDataFrame(df_annot_muscle,
                                           keep.extra.columns = T,
                                           seqnames.field="seqnames")
gr_annot_liver = makeGRangesFromDataFrame(df_annot_liver,
                                          keep.extra.columns = T,
                                          seqnames.field="seqnames")

gr_full_muscle = join_overlap_full(gr_annot_muscle,gr_QTLdb_filtered)
gr_full_filtered_muscle = subsetByOverlaps(subsetByOverlaps(gr_full_muscle, gr_annot_muscle), gr_QTLdb_filtered)

gr_left_muscle = subsetByOverlaps(gr_full_muscle, gr_annot_muscle)
df_left_muscle =as.data.frame(gr_left_muscle)
df_left_muscle$cat_trait.y[is.na(df_left_muscle$cat_trait.y)]="None"
contigency_table = table(df_left_muscle$cat_trait.y, df_left_muscle$Annotation.x)

df_fisher_muscle = data.frame(Category = rownames(contigency_table),
                              pval = rep(1, 6),
                              OR = rep(0,6),
                              IC_min = rep(0,6),
                              IC_max = rep(0,6))
df_fisher_muscle$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$p.value})
df_fisher_muscle$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$estimate})
df_fisher_muscle$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[1]})
df_fisher_muscle$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[2]})


df_QTL_muscle = as.data.frame(gr_full_filtered_muscle)
df_QTL_muscle$Freq = rep(1, nrow(df_QTL_muscle))
freq_QTL_annot_muscle = aggregate(df_QTL_muscle$Freq,
                                  list(df_QTL_muscle$Annotation.x, df_QTL_muscle$cat_trait.y),
                                  sum)
colnames(freq_QTL_annot_muscle) = c("Annotated","Trait category","Total")
freq_QTL_annot_muscle$Freq = freq_QTL_annot_muscle$Total / rep(table(df_annot_muscle$Annotation),5)


gr_full_liver = join_overlap_full(gr_annot_liver,gr_QTLdb_filtered)
gr_full_filtered_liver = subsetByOverlaps(subsetByOverlaps(gr_full_liver, gr_annot_liver), gr_QTLdb_filtered)

gr_left_liver = subsetByOverlaps(gr_full_liver, gr_annot_liver)
df_left_liver =as.data.frame(gr_left_liver)
df_left_liver$cat_trait.y[is.na(df_left_liver$cat_trait.y)]="None"
contigency_table = table(df_left_liver$cat_trait.y, df_left_liver$Annotation.x)

df_fisher_liver = data.frame(Category = rownames(contigency_table),
                             pval = rep(1, 6),
                             OR = rep(0,6),
                             IC_min = rep(0,6),
                             IC_max = rep(0,6))
df_fisher_liver$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$p.value})
df_fisher_liver$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$estimate})
df_fisher_liver$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[1]})
df_fisher_liver$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),1:2])$conf.int[2]})

df_fisher_REG = data.frame(Category = rownames(contigency_table),
                           pval = rep(1, 6),
                           OR = rep(0,6),
                           IC_min = rep(0,6),
                           IC_max = rep(0,6))
df_fisher_REG$pval = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$p.value})
df_fisher_REG$OR = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$estimate})
df_fisher_REG$IC_min = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$conf.int[1]})
df_fisher_REG$IC_max = sapply(1:6, function(x) {fisher.test(contigency_table[c(x,4),3:2])$conf.int[2]})

df_fisher_all = rbind(cbind(df_fisher_muscle, Tissue = "Muscle"),
                      cbind(df_fisher_liver, Tissue = "Liver"),
                      cbind(df_fisher_REG, Tissue = "Not tissue specific"))

df_fisher_all$Threshold = "False"
df_fisher_all$Threshold[which(df_fisher_all$pval<0.05)]="True"

df_fisher_all %>% dplyr::select(-Threshold) %>%
    dplyr::select(Tissue, everything()) %>%
    mutate(pval=round(pval, 3), OR = round(OR, 3), IC_min = round(IC_min, 3),
           IC_max= round(IC_max, 3)) %>% View



## Obsolete --------------------------------------------------------------------

### Another --------------------------------------------------------------------

df_cor_Vbeta_tissue = data.frame(Scenario = rep(NA, 30),
                                 Cor_Vbeta = rep(NA, 30),
                                 Cor_Vbeta_sum = rep(NA, 30))
compt = 1
for (i in setdiff(unique(Vbeta_annot$gene),"LEPR")){
    temp = filter(Vbeta_annot, gene== i)
    df_cor_Vbeta_tissue$Scenario[compt] = paste0(i, "_DU",sep="")
    df_cor_Vbeta_tissue$Scenario[compt+1] = paste0(i, "_LW",sep="")
    df_cor_Vbeta_tissue$Scenario[compt+2] = paste0(i, "_LD",sep="")
    df_cor_Vbeta_tissue$Cor_Vbeta_sum[compt] = cor(as.numeric(filter(temp, learning=="DU" & tissue =="liver")$Sum_vi),
                                                   as.numeric(filter(temp, learning=="DU" & tissue == "muscle")$Sum_vi))
    df_cor_Vbeta_tissue$Cor_Vbeta_sum[compt+1] = cor(as.numeric(filter(temp, learning=="LW" & tissue =="liver")$Sum_vi),
                                                     as.numeric(filter(temp, learning=="LW" & tissue == "muscle")$Sum_vi))
    df_cor_Vbeta_tissue$Cor_Vbeta_sum[compt+2] = cor(as.numeric(filter(temp, learning=="LD" & tissue =="liver")$Sum_vi),
                                                     as.numeric(filter(temp, learning=="LD" & tissue == "muscle")$Sum_vi))
    df_cor_Vbeta_tissue$Cor_Vbeta[compt] = cor(as.numeric(filter(temp, learning=="DU" & tissue =="liver")$Vbeta),
                                               as.numeric(filter(temp, learning=="DU" & tissue == "muscle")$Vbeta))
    df_cor_Vbeta_tissue$Cor_Vbeta[compt+1] = cor(as.numeric(filter(temp, learning=="LW" & tissue =="liver")$Vbeta),
                                                 as.numeric(filter(temp, learning=="LW" & tissue == "muscle")$Vbeta))
    df_cor_Vbeta_tissue$Cor_Vbeta[compt+2] = cor(as.numeric(filter(temp, learning=="LD" & tissue =="liver")$Vbeta),
                                                 as.numeric(filter(temp, learning=="LD" & tissue == "muscle")$Vbeta))
    compt = compt + 3
}

df_cor_Vbeta_tissue$Gene = rep(setdiff(unique(Vbeta_annot$gene),"LEPR"),each=3)
df_cor_Vbeta_tissue$Breed = rep(c("DU","LW","LD"))

ggplot(df_cor_Vbeta_tissue,aes(x=Cor_Vbeta, color = Breed,fill=Breed))+
    geom_density(alpha=0.3)+
    theme_bw()+
    xlab("Correlation between individual Vi in liver and muscle")+
    xlim(0,1)

ggplot(df_cor_Vbeta_tissue,aes(x=Cor_Vbeta_sum, color = Breed,fill=Breed))+
    geom_density(alpha=0.3)+
    theme_bw()+
    xlab("Correlation between window-based Vi in liver and muscle")+
    xlim(0,1)

## Option 2
# fig5b_option2 <- ggplot(results %>% filter(learning != validation,
#                           annotation %in% c("none", "vep+atac+methyl"))) +
#     geom_point(aes(x=annotation, y=cor_spearman, color=learning)) +
#     geom_line(aes(x=annotation, y=cor_spearman, color=learning,
#                   group=paste0(learning, ".", validation))) +
#     facet_grid(tissue ~ Gene.Name) +
#     geom_hline(yintercept = 0, lty = 2) +
#     theme_bw() +
#     ylim(c(-.3,1)) +
#     theme(axis.text.x = element_text(angle = 45, hjust=1))

## Option 3: BayesRCO - BayesR
# results_summary <- results %>%
#     filter(annotation %in% c("none", "vep+atac+methyl"),
#            learning != validation) %>%
#     dplyr::select(-cor_pearson, -Chromosome, -EnsemblID) %>%
#     group_by(learning, validation, tissue, Gene.Name) %>%
#     summarize(spearman_diff = diff(cor_spearman))
# fig5b_option3 <- ggplot(results_summary) +
#     geom_point(aes(x=Gene.Name, y=spearman_diff,
#                    color=paste0(learning, ".", validation))) +
#     facet_grid(~tissue) +
#     geom_hline(yintercept = 0, lty = 2) +
#     theme_bw() +
#     xlab("") + ylab("BayesRCO-BayesR validation Spearman correlation") +
#     theme(axis.text.x = element_text(angle = 45, hjust=1))

### Per-chr genomic PCAs by annotation  ----------------------------------------
# eigenvec <- vector("list", length=8)
# names(eigenvec) <- paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18))
# eigenval <- vector("list", length=8)
# names(eigenval) <- paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18))
# for(annot in c("general", "atacseq.methylation3", "atacseq", "methylation3", "vep")) {
#
#     for(chr in names(eigenvec)) {
#         eigenvec[[chr]] <- read.table(paste0("pca/WGS_300_GS_filter_",
#                                              chr, "_", annot, ".eigenvec"))
#         eigenval[[chr]] <- read.table(paste0("pca/WGS_300_GS_filter_",
#                                              chr, "_", annot, ".eigenval"))
#     }
#
#     eigenvec_df <- bind_rows(eigenvec, .id="chr")
#     colnames(eigenvec_df)[c(2:5)] <- c("breed", "ID", "PCA1", "PCA2")
#     eigenvec_df$chr <- factor(eigenvec_df$chr,
#                               levels = paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18)))
#
#     perc1 <- lapply(eigenval, function(x) {
#         data.frame(PC1 = paste0("PC1: ",round(x$V1[1] / sum(x$V1)*100,1), "%"))}) %>%
#         bind_rows(., .id="chr")
#     perc2 <- lapply(eigenval, function(x) {
#         data.frame(PC2 = paste0("PC2: ",round(x$V1[2] / sum(x$V1)*100,1), "%"))}) %>%
#         bind_rows(., .id="chr")
#     perc <- full_join(perc1, perc2, by="chr")
#     perc$perc <- paste0(perc$PC1, "\n", perc$PC2)
#     perc$chr <- factor(perc$chr,
#                        levels = paste0("chr", c(1, 2, 5, 6, 7, 10, 14, 18)))
#
#     annot_title <- ifelse(annot == "general", "All annotations",
#                           ifelse(annot == "atacseq", "ATAC-seq regions",
#                                  ifelse(annot == "methylation3", "Methylation regions",
#                                         ifelse(annot == "vep", "REG+NSC",
#                                                "ATAC-seq and methylation regions"))))
#     g <- ggplot(eigenvec_df) +
#         geom_point(aes(x=PCA1, y=PCA2, col=breed), alpha = 0.3) +
#         geom_hline(aes(yintercept = 0), lty=2, col="grey30") +
#         geom_vline(aes(xintercept = 0), lty=2, col="grey30") +
#         geom_text(data=perc, aes(x=Inf, y=-Inf, label=perc), hjust = 1, vjust = -.1) +
#         facet_wrap(~chr) +
#         theme_bw() +
#         scale_color_discrete(name="") +
#         theme(legend.position="bottom") +
#         ggtitle(annot_title)
#
#     print(g)
#
# }
