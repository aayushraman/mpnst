### Heatmap for Single ChIP-seq Mark ##
rm(list = ls())
library(AnnotationDbi)
library(DiffBind)
library(ComplexHeatmap)
library(dplyr)
library(genefilter)
library(ggplot2)
library(GenomicRanges)
library("org.Hs.eg.db")
library(pheatmap)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


setwd("/Volumes/My_Passport/MPNST_cell_lines_compiled_Chip_seq_data/04aln_downsample/")
H3K27me3 <- dba(sampleSheet = "H3K27me3.csv")
H3K27me3_reads_minus <- dba.count(H3K27me3, fragmentSize = 200, score = DBA_SCORE_READS_MINUS)
save(H3K27me3_reads_minus, file = "H3K27me3_reads_minus.RData")
load("/Volumes/My_Passport/MPNST_cell_lines_compiled_Chip_seq_data/04aln_downsample/H3K27me3_reads_minus.RData")


## from dba count object
reads.RPKM <- dba.peakset(H3K27me3_reads_minus, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
name <- "RPKM"
rownames(reads.RPKM) <- paste(reads.RPKM$CHR,reads.RPKM$START,reads.RPKM$END,sep = ":")
reads.RPKM.mat <- reads.RPKM[, 4:ncol(reads.RPKM)] %>% as.matrix()
colnames(reads.RPKM.mat) <- as.vector(H3K27me3_reads_minus$class["Tissue",])
genotypes <- factor(c(rep("PRC2 loss",5),rep("PRC2 WT", 4), rep("NSC")), levels = c("PRC2 loss", "PRC2 WT", "NSC"))
sample.lab <- as.vector(H3K27me3_reads_minus$class["Tissue",])

## All promoters
reads.RPKM.mat.promoter <- reads.RPKM.mat

## Top regions
ind <- which(reads.RPKM.mat.promoter <  0)
reads.RPKM.mat.promoter[ind] <- 0
hmcols <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
# all.peaks <- varFilter(eset = reads.RPKM.mat.genebody, var.cutoff = 0.75) ## Top 10000 peaks
all.peaks <- CancerSubtypes::FSbyMAD(Data = reads.RPKM.mat.promoter, cut.type = "topk", value = 10000)
print(dim(all.peaks))
plot(H3K27me3_reads_minus)


## MDS Plot
mds.plot <- cmdscale(d = dist(t(log2(all.peaks+1))), k = 2)
ggplot(data = data.frame(mds.plot), aes(mds.plot[,1],mds.plot[,2], color = genotypes)) +
  xlab("Comp 1") + ylab("Comp 2") + geom_point(size = 4) + 
  ggrepel::geom_text_repel(aes(label = sample.lab)) + 
  theme(legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        plot.title = element_text(size = 0, face ="bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

## Correlation plot
annot1 <- data.frame(genotypes, row.names = sample.lab)
annot1$genotypes <- factor(annot1$genotypes, levels = c("PRC2 WT", "PRC2 loss", "NSC"))
pheatmap(cor(all.peaks), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE,
         show_colnames = TRUE, fontsize_row = 10, legend = TRUE, filename = NA,
         fontsize_col = 10, scale = "none", fontface="bold",
         clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
         annotation = annot1)
## heatmap
annot1 <- data.frame(genotypes, row.names = sample.lab)
annot1$genotypes <- factor(annot1$genotypes, levels = c("NSC", "PRC2 WT", "PRC2 loss"))
pheatmap(log2(all.peaks+1), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
         show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA,
         fontsize_col = 10, scale = "row", fontface="bold",
         clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
         annotation = annot1)
dba.plotHeatmap(H3K27me3_reads_minus, correlations=FALSE)

## GREAT
file <- "H3K27me3"
file <- paste("~/Desktop/GREAT_",file,"All_diff.feat.selected.txt",sep = "") ## change the directory if you want
bedfile <- data.frame()
for(j in 1:nrow(all.peaks)){
  bedfile[j, 1] = noquote(strsplit(rownames(all.peaks)[j], split = ":|-")[[1]])[1]
  bedfile[j, 2] = noquote(strsplit(rownames(all.peaks)[j], split = ":|-")[[1]])[2]
  bedfile[j, 3] = noquote(strsplit(rownames(all.peaks)[j], split = ":|-")[[1]])[3]
}
colnames(bedfile) = c("chrom", "chromStart", "chromEnd")
write.table(x = bedfile, file = file, sep="\t", quote=F, row.names = FALSE, col.names = FALSE)


## background file
background <- data.frame(regions = rownames(reads.RPKM.mat))
background$chr <- gsub(pattern = "(.*):(.*):(.*)", replacement = "\\1", x = background$regions)
background$start <- gsub(pattern = "(.*):(.*):(.*)", replacement = "\\2", x = background$regions)
background$end <- gsub(pattern = "(.*):(.*):(.*)", replacement = "\\3", x = background$regions)
write.table(x = background[,2:4], file = "~/Desktop/All_Background.txt", sep="\t", quote=F, 
            row.names = FALSE, col.names = FALSE)

### Supervised Heatmaps from Differential Analysis ###

## Contrast, analyze and report from DiffBind ##
H3K27me3_contrast <- dba.contrast(H3K27me3_reads_minus, 
                                  H3K27me3_reads_minus$masks$PRC2loss,
                                  H3K27me3_reads_minus$masks$PRC2WT,
                                  "PRC2_loss", "PRC2_WT") ## CHANGE
H3K27me3_analyze <- dba.analyze(H3K27me3_contrast)
dba.report(H3K27me3_analyze)
plot(H3K27me3_analyze, contrast=1)
plot(H3K27me3_analyze, contrast=1, correlations=FALSE)
peaks <- data.frame(dba.report(H3K27me3_analyze))
peaks_up <- peaks[which(peaks$Fold > 0 & peaks$FDR < 0.05),]
peaks_down <- peaks[which(peaks$Fold < 0 & peaks$FDR < 0.05),]
write.table(x = peaks_up, file = "~/Desktop/H3K27me3_Peaks_H3K27me3_loss_up.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = TRUE)
write.table(x = peaks_up[,c(1:3)], file = "~/Desktop/H3K27me3_Peaks_H3K27me3_loss_up-GREAT.bed", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
dim(peaks_up)
write.table(x = peaks_down, file = "~/Desktop/H3K27me3_Peaks_H3K27me3_loss_down.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = TRUE)
write.table(x = peaks_down[,c(1:3)], file = "~/Desktop/H3K27me3_Peaks_H3K27me3_loss_down-GREAT.bed", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
dim(peaks_down)

## Supervised Heatmaps
H3K27me3_analyze <- dba.analyze(H3K27me3_contrast)
degs <- data.frame(dba.report(H3K27me3_analyze))
rownames(degs) <- paste(degs$seqnames,":",degs$start,"-",degs$end,sep = "")
plot(H3K27me3_analyze, contrast=1)
mat <- plot(H3K27me3_analyze, contrast=1, correlations=FALSE)
mat <- data.frame(mat)
rownames(mat) <- paste(mat$seqnames,":",mat$start,"-",mat$end,sep = "")
col.names <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
mat <- inner_join(degs %>% tibble::rownames_to_column(var = "region"),
                  mat[,col.names] %>% tibble::rownames_to_column(var = "region"), by = "region")
annot1 <- data.frame(genotypes, row.names = sample.lab)
annot1$genotypes <- factor(annot1$genotypes, levels = c("PRC2 WT", "NSC", "PRC2 Neg"))
colnames(mat)[13:19] <- sample.lab
rownames(mat) <- mat$region
pheatmap(mat[,13:19], cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
         show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA,
         fontsize_col = 10, scale = "row", fontface="bold", clustering_distance_rows = "correlation", 
         clustering_distance_cols = "euclidean", annotation = annot1)
dba.plotHeatmap(H3K27me3_analyze, correlations=FALSE, contrast = 1)

