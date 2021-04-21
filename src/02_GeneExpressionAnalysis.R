#############################################################################
# @Author: Ayush T. Raman
# Rai Lab, MD Anderson Cancer Center
# Date: 10th Dec. 2016
# Updated: v2: 30th Oct 2017
#
# Program is used for:
# 1. Diff Expression Analysis for RNA-Seq Data (Final Analysis)
# 2. similar code was used for the other cell line datasets
#
# Comments
# 1. Comparison of MPNST vs Normals only
#############################################################################

## Clearing Workspace
rm(list = ls())

## libraries used
library(DESeq2)
library(dplyr)
library(genefilter)
library(ggplot2)
library(ggrepel)
library(pcaExplorer)
library(pheatmap)
library(rafalib)
library(RColorBrewer)
library(tibble)
library(Vennerable)
mypar()

## source all the function
source("combiningFilesinDirectory.R")
source("DESeqCalculation.R")
source("plots2Plot.R")
source("tumor.comp_plots.R")

############################
#
## Data Analysis
#
############################

## files
filePath <- "counts/"
file_list <- list.files(path = filePath, pattern = "(.*)ReadsPerGene.out.tab")
counts.table <- combiningFilesinDirectory(filePath = filePath, file_list = file_list, n = 2)
counts.table <- counts.table[-c(1:3),]

## All samples -- prepare the data matrix
mpnst <- counts.table %>% select(matches("(T3|T6|T9|32T|7BT|T26|T27)")) %>% colnames()
mpnst <- mpnst[c(7:13)]
normal <- counts.table %>% select(matches("Sample_(.*)N([0-9]*)ReadsPerGene.out.tab")) %>% colnames() ## 13 samples
normal <- normal[which(!normal %in% c("Sample_14NReadsPerGene.out.tab", "Sample_NF5NReadsPerGene.out.tab",
                                      "Sample_25NReadsPerGene.out.tab", "Sample_N1ReadsPerGene.out.tab",
                                      "Sample_31NReadsPerGene.out.tab", "Sample_NF4NReadsPerGene.out.tab"))]
counts.neuro.mpnst <- counts.table[,c(normal, mpnst)]
counts.neuro.mpnst <- counts.neuro.mpnst[,c(1:7,11:17)]
dim(counts.neuro.mpnst)
colnames(counts.neuro.mpnst)
genotypes <- factor(c(rep("Normal", 7), rep("MPNST",7)), levels = c("MPNST", "Normal"))
colnames(counts.neuro.mpnst) <- gsub(pattern = "Sample_(.*)ReadsPerGene.out.tab",
                                     replacement = "\\1", x = colnames(counts.neuro.mpnst))
colData <- data.frame(colnames(counts.neuro.mpnst), genotypes)
rownames(colData) <- as.character(colData[,1])
head(counts.neuro.mpnst)
head(colData)

## DESeq2 results MPNST vs Normal
dim(counts.neuro.mpnst)
dds <- DESeqDataSetFromMatrix(counts.neuro.mpnst, colData, design = ~ genotypes)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
idx <- which(rowMeans(dat) > 30)
dat <- dat[idx,]
dds <- dds[idx,]
print(dim(dat))
print(dim(dds))

## size factor
sizeFactorRatio <- sizeFactors(dds)
readCounts <- colSums(counts(dds))
qplot(x = sizeFactors(dds), y = colSums(counts(dds))/1e6) + geom_point(aes(colour = dds$genotypes), size = 5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30") + 
  geom_text_repel(aes(label = colnames(assay(dds))), size = 5, fontface = "bold", color = "black", 
                  box.padding = unit(0.35, "lines"), point.padding = unit(0.2, "lines")) + 
  xlab("Size Factor") + ylab("Number of Aligned Reads (in million)") +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), legend.position = c(.85, .1),
        legend.title = element_blank(), legend.text = element_text(size = 14, face = "bold", color = "black"))

## PCA Plot
rld.dds <- rlog(dds, blind = FALSE)
print(plotPCA(rld.dds, intgroup = c("genotypes")))
p <- plotPCA(rld.dds, intgroup = c("genotypes")) 
p + geom_text_repel(aes(label = colnames(assay(rld.dds))), size = 5, 
                    fontface = "bold", color = "black", box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.2, "lines")) + coord_fixed() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        legend.title = element_blank(), legend.position = c(.85,.1),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
p1 <- PCAplot(data = assay(rld.dds), genotypes = colData$genotypes, conditions = colData$genotypes, title = "")
p1$plot + geom_text_repel(aes(label = colnames(assay(rld.dds))), size = 5, 
                          fontface = "bold", color = "black", box.padding = unit(0.35, "lines"),
                          point.padding = unit(0.2, "lines")) + coord_fixed() + theme_bw(base_size = 16) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        legend.title = element_text(size = 18, face = "bold", color = "black"),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

## housekeeping genes
genes <- c("GAPDH", "ACTB", "ATF2", "PPP2R5B", "AKT1", "EZH2", "IGF2BP1", "MYL2")
plotBox_genes(counts.table = log2(dat+1), list.Gname = genes, genotypes = genotypes, geneNames = genes)

## PCA -- using the rld object
pcs <- prcomp(t(assay(rld.dds)), center = TRUE)
percentVar <- round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2)

###########################
#
## DESeq Analysis
#
###########################

dds <- DESeq(dds, betaPrior = TRUE)
resultsNames(dds)
colnames(colData)[1] <- "Samples"
res.dds <- results(dds, contrast = c("genotypes","MPNST", "Normal"))
ind.mpnst <- c(1:14)
mpnst.normal <- plots.mpnst(dds = dds, res.dds = res.dds, cols = ind.mpnst, 
                            fc = 1.5, lev1 = c(1,2), geno.cols = c("#7570b3", "#1b9e77"), 
                            tr.fal = FALSE, genotypes = genotypes, rld = rld.dds)
mpnst.normal$plot.all
mpnst.normal$plot.up
mpnst.normal$plot.down
mpnst.normal$plot.Wnt.Cat.withoutlabels + coord_cartesian(ylim = c(0,35)) + xlab("Log2 Fold Change") + ylab("-Log10(P-value)")
mpnst.normal$plot.Wnt.Cat + coord_cartesian(xlim = c(-7.5,10), ylim = c(0,25)) + 
    xlab("Log2 Fold-Change") + ylab("-Log10 (P-value)") + theme(legend.position = "none") 
save(mpnst.normal, file = "MPNST.Tumor.RData")

##################################
#
## Plots -- Gene Expression
#
##################################


load("MPNST.Tumor.RData")
library(reshape2)
library(ggpubr)
neural.genes <- c("SOX9", "MSX1", "ZIC1", "ZIC2", "ZIC3",
                  "DLX5", "DLX6", "SNAI2", "ID1", "PAX7", "SOX5", "PAX3")
schwann.genes <- c("MBP", "S100B", "MPZ", "L1CAM", "NCAM2", "GAP43", "SOX2")
neural.genes <- schwann.genes
mat <- mpnst.normal$result$results[mpnst.normal$result$results$gene %in% 
                                       neural.genes,c(7:20,22)]
mat <- mat[order(mat$gene), ]
mat <- melt(mat)
mat$genotypes <- rep(genotypes, each = length(neural.genes))
mat$value <- log2(mat$value + 1)
mat$genotypes <- relevel(mat$genotypes, ref = "Normal")
p1 <- ggboxplot(mat, "gene", "value", fill = "genotypes",
                palette = c("#00cccc","#ff4d4d"), notch = TRUE) + theme_grey() +
        facet_grid(.~ gene,  space = "free", scale = "free") + ylab("Log2 Normalized Counts") + 
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
              axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
              axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.title.x =  element_blank(), axis.line.x = element_blank(),
              legend.title = element_blank(), 
              legend.text = element_text(size = 18, face = "bold", color = "black"),
              strip.text = element_text(size = 14, colour = "black", face = "bold"))
filename <- "/Users/araman/Desktop/MPNST/results_Tumor-CellLines-final/SchwannCell.png"
ggsave(filename = filename, plot = p1, scale = 1, width = 8, height = 6,
       units = c("in"), dpi = 600, limitsize = TRUE)

## Plots for EZH2 and CDH3
genes <- c("CDH3","EZH2")
dat1 <- mpnst.normal$result$results[mpnst.normal$result$results$gene %in% 
                                       genes,]
for(i in 1:length(genes)[1]){
    mat <- data.frame(labels = factor(colnames(dat1[,c(7:20)]), 
                                      levels = colnames(dat1[,c(7:20)])),
                      Normalized.Counts = as.vector(as.matrix(dat1[i,c(7:20)])), 
                      Genotypes = genotypes)
    mat$Genotypes <- relevel(mat$Genotypes, ref = "Normal")
    pval <- dat1[i,"padj"]
    logFC <- dat1$log2FoldChange[i]
    pval.txt <- paste("FDR = ", format(x = pval, digits = 4), 
                      " & LogFC = ", format(x = logFC, digits = 4),sep = "")
    message(pval.txt)
    mat$Normalized.Counts <- log2(mat$Normalized.Counts + 1)
    lab.txt <- max(mat$Normalized.Counts) + 1.5
    p1 <- ggviolin(mat, x = "Genotypes", y = "Normalized.Counts", fill = "Genotypes",
                   palette = c("#00AFBB","#FC4E07"), add = "boxplot", 
                   add.params = list(fill = "white")) +
        # stat_compare_means(data = mat[-c(4,10),], comparisons = list(c("MPNST", "Normal")), 
        #                    method = "t.test", label = "p.signif", paired = TRUE) + 
        xlab("") + ylab("Log2 Normalized Counts") + ylim(c(-3,14)) +
        annotation_custom(grob = text_grob(pval.txt, face = "bold", hjust = 0.5, vjust = 0.5),
                          xmin = 1, xmax = 2, ymin = lab.txt-3.25, ymax = lab.txt+4) +
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.title.y= element_text(size = 16, colour = "black", face = "bold"),
              axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
              axis.ticks.x = element_blank(), axis.title.x =  element_blank(), 
              axis.line.x = element_blank(), axis.text.x =  element_blank(),
              legend.title = element_blank(), 
              legend.text = element_text(size = 14, face = "bold"))
    filename <- paste("MPNST/results_Tumor-CellLines-final/",
                      genes[i],".png",sep = "")
    ggsave(filename = filename, plot = p1, scale = 1, width = 6, height = 4.5,
           units = c("in"), dpi = 500, limitsize = TRUE)
}


## Wnt.Cat
WNT.Cat <- read.table("WNT-Cat_Patients.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
sum(WNT.Cat$V1 %in% rownames(assay(dds)))
g1 <- WNT.Cat$V1[which(!WNT.Cat$V1 %in% rownames(assay(dds)))]
g1 %in% rownames(counts.table)
rm(g1)
Wnt.Cat.tumors <- WNT.Cat$V1[WNT.Cat$V1 %in% rownames(mpnst.normal$result$up.reg)]

## Comparison with ST88/S462 and Tumor Sample
load(file = "MPNST/ST88/results.ST88.RData")
results.ST88$results <- results.ST88$results[!is.na(results.ST88$results$padj),]
plot1 <- volcano.plot(results = results.ST88$results, fc = 1, pval = 0.05)
plot1$plot.Wnt.Cat.withoutlabels + coord_cartesian(ylim = c(0,80)) + xlab("Log2 Fold Change") + ylab("-Log10(P-value)")
plot1$plot.Wnt.Cat + coord_cartesian(ylim = c(0,80)) + xlab("Log2 Fold Change") + ylab("-Log10(P-value)")

load(file = "MPNST/S462/results.S462.RData")
results.S462$results <- results.S462$results[!is.na(results.S462$results$padj),]
plot2 <- volcano.plot(results = results.S462$results, fc = 1.5, pval = 0.05)
plot2$plot.Wnt.Cat.withoutlabels + xlab("Log2 Fold Change") + ylab("-Log10(P-value)")
plot2$plot.Wnt.Cat + xlab("Log2 Fold Change") + ylab("-Log10(P-value)")

## Venn Diagram and Common Genes
Wnt.Cat.ST88 <- WNT.Cat$V1[WNT.Cat$V1 %in% rownames(results.ST88$up.reg)]
Wnt.Cat.S462 <- WNT.Cat$V1[WNT.Cat$V1 %in% rownames(results.S462$up.reg)]
w1 <- Venn(Sets=list(Tumors = WNT.Cat$V1, ST88 = Wnt.Cat.ST88, S462 = Wnt.Cat.S462))
plot(w1, doWeights=F)
common.genes <- WNT.Cat[Wnt.Cat.tumors %in% Wnt.Cat.ST88]
common.genes[common.genes %in% Wnt.Cat.S462]

## Heatmap -- ST88, S462 and NSC
mat1 <- results.ST88$results[which(rownames(results.ST88$results) %in% WNT.Cat$V1),c(7:12)]
mat1$genes <- rownames(mat1)
mat2 <- results.S462$results[which(rownames(results.S462$results) %in% WNT.Cat$V1),c(7:12)]
mat2$genes <- rownames(mat2)
mat <- inner_join(mat1, mat2, by = "genes")
head(mat)
rownames(mat) <- mat$genes
colnames(mat) <- gsub(pattern = "norm.counts.", replacement = "", x = colnames(mat))
mat <- mat[,c("Control_NSC_1.x","Control_NSC_2.x","Control_NSC_3.x",
              "ST88_1","ST88_2","ST88_3","S462_1","S462_2","S462_3")]
genotypes <- factor(c(rep("NSC",3),rep("ST88",3),rep("S462",3)),levels = c("NSC","ST88","S462"))
annot1 <- data.frame(genotypes, row.names = colnames(mat))
pheatmap(mat = log2(mat+1), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE,
         show_colnames = FALSE, legend = TRUE, filename = NA, fontsize = 10, scale = "row", 
         fontface="bold", annotation = annot1)


## X-Y plots for proof of concept
# source('~/Desktop/MPNST/src/comp.logFC.R')
# dat2 <- results.ST88$results
# dat2$gene <- rownames(results.ST88$results)
# comp.logFC(dat1 = mpnst.normal$result$results, dat2 = dat2, sample1 = "Tumors", sample2 = "ST88")

## Heatmap on the types of genes
mat <- mpnst.normal$result$results
genes <- c("TWIST1", "SNAI2", "CDH2", "FN1", "COL1A2", "COL5A2", "FOXC2", "GSC", "MMP9", "STEAP1",
           "TIMP1", "TMEM132A", "WNT5A", "BMP7", "GATA4", "BMP1", "SOX9", "ITGA4", "CHD7",
           "MBP", "S100B", "PROM1", "SOX2")
mat <- mat[match(genes, mat$gene), c(7:20)]
genotypes <- factor(c(rep("Normal", 7), rep("MPNST",7)), levels = c("Normal","MPNST"))
annot1 <- data.frame(genotypes, row.names = colnames(mat))
pheatmap(mat = mat, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE,
         show_colnames = FALSE, legend = TRUE, filename = NA, fontsize = 10, scale = "row", 
         fontface="bold", annotation = annot1)

## Using ComplexHeatmap
df1 <- data.frame(Samples = genotypes)
geno.cols <- c("#00e5e5", "#e50000")
geno.cols.assigned <- setNames(geno.cols, unique(levels(df1$Samples)))
annot1 <- HeatmapAnnotation(df = df1, col = list(Samples = geno.cols.assigned), show_legend = TRUE, 
                            annotation_legend_param = list(title = "Genotypes", width = unit(1, "cm"), 
                                                           title_gp = gpar(fontsize = 14,fontface = "bold"), 
                                                           labels_gp = gpar(fontsize = 14, fontface = "bold")))
rld.dds <- rld.dds[match(genes, rownames(rld.dds)),]
mat.scaled <- t(scale(t(assay(rld.dds)),center=TRUE,scale=FALSE))
Heatmap(mat.scaled, show_row_names = TRUE, show_column_names = FALSE, row_dend_reorder = TRUE, 
            column_dend_reorder = TRUE, clustering_distance_columns = "euclidean", 
            cluster_rows = FALSE, row_names_gp = gpar(fontsize = 14, fontface = "bold"),
            column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
            clustering_method_columns = "complete", top_annotation = annot1,
            heatmap_legend_param = list(title = "Scaled Counts",
                                        title_gp = gpar(fontsize = 14, fontface = "bold"),
                                        labels_gp = gpar(fontsize = 14, fontface = "bold")))
