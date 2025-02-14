#library(ggplot2)
library(dplyr, lib="/home/araul/R/3.0/library/")
library(tidyr, lib="/home/araul/R/3.0/library/")
library(GenomicRanges, lib="/home/araul/R/3.0/library/")

#BIMPATH=$2         #bim file with coordinates (same order than param)
#QTLDBPATH=$3       #QTLdb file
#WRITEPATH=$4        #output
#extend_bp=$5       #windows around QTL, 0 by default, will be applied both start and end of the QTL

# inputs
QTLDBPATH="/travail/araul/GS_data/QTLdb_pigSS11.bed"    #QTLdb file
extend_bp=0         #windows around QTL, 0 by default, will be applied both start and end of the QTL

QTLdb <- read.table(QTLDBPATH,
                                        skip = 18,
                                        fill = T,
                                        row.names = NULL,
                                        header = F,
                                        sep="\t")
QTLdb = QTLdb[,-c(5:12)]
colnames(QTLdb) = c("chromosome","start","end","trait")
QTLdb$trait = sub("\\(.*", "", QTLdb$trait)
QTLdb$chromosome = sub("....","", QTLdb$chromosome)
QTLdb$start = QTLdb$start-extend_bp
QTLdb$end = QTLdb$end+extend_bp
gr_QTLdb = makeGRangesFromDataFrame(QTLdb, keep.extra.columns = T)

# annex function
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
      mcols(uniqex) <- mcols(uniqex)[,c((length(cns.x)+1):(length(cns.x)+length(cns.y)),1:(length(cns.x)))]
      mcols(uniqey)[, (length(cns.y)+1):(length(cns.x)+length(cns.y))] <- NA
      colnames(mcols(uniqey)) <- c(cns.y,cns.x)
      res <- c(z, uniqex,uniqey)
      mcols(res) <- mcols(res) %>% data.frame() %>% mutate_all(as.character())
      return(res)
  }

gr_full_filtered_list_none <- list()
index_none <- 1
gr_full_filtered_list_annot <- list()
index_annot <- 1

for(PARAMPATH in dir("/travail/araul/GS_data/acrossbreed_results/")) {  #param file
      annot <- unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[5]))
      cat(PARAMPATH, "\n")
      chr <- unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[1]))
      BIMPATH <- paste0("/travail/araul/GS_data/WGS_300_GS_filter_", chr, ".bim")

      # read files
      paramfile <- read.table(paste0("/travail/araul/GS_data/acrossbreed_results/", PARAMPATH), header=T)
      bimfile = read.table(BIMPATH)
      bimfile_v2 = bimfile[,c(1,2,4,4)]
      colnames(bimfile_v2) = c("chromosome","id","start","end")
      gr_bim = makeGRangesFromDataFrame(bimfile_v2, keep.extra.columns = T)
      param_bim=cbind(paramfile,bimfile_v2)
      gr_param_bim = makeGRangesFromDataFrame(param_bim, keep.extra.columns = T)

      # overlap with granges
      gr_full = join_overlap_full(gr_param_bim,gr_QTLdb)
      gr_full_filtered = subsetByOverlaps(subsetByOverlaps(gr_full, gr_param_bim), gr_QTLdb)
      gr_full_filtered_df <- as.data.frame(gr_full_filtered)
      gr_full_filtered_df$gene <- unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[2]))
      gr_full_filtered_df$tissue <- unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[3]))
      gr_full_filtered_df$annotation <- unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[5]))
      gr_full_filtered_df$learning <- unlist(strsplit(unlist(lapply(strsplit(PARAMPATH, split="_", fixed=TRUE), function(x) x[6])), split=".param", fixed=TRUE))
      if(annot == "none") {
            gr_full_filtered_list_none[[index_none]] <- gr_full_filtered_df
                index_none <- index_none + 1
          } else {
                gr_full_filtered_list_annot[[index_annot]] <- gr_full_filtered_df
                index_annot <- index_annot + 1
              }
  }

annot_results <- do.call("rbind", gr_full_filtered_list_annot)
none_results <- do.call("rbind", gr_full_filtered_list_none)

#write param files filtered by co-localization with QTL
write.table(annot_results,
                        file="/travail/araul/GS_data/coloc_param-QTL_annot.txt",
                        sep="\t",
                        row.names=F,
                        col.names=T,
                        quote=F)

write.table(none_results,
                        file="/travail/araul/GS_data/coloc_param-QTL_none.txt",
                        sep="\t",
                        row.names=F,
                        col.names=T,
                        quote=F)
