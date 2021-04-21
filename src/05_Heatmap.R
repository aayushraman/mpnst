############################################################
# @Author: Ayush T. Raman
# Rai Lab, MD Anderson Cancer Center
# Date: Feb 26th, 2019
#
# Program is used for:
# 1. Common Genes across Celllines and Tumor Datasets

## Clearing Workspace
rm(list = ls())

## libraries used
library(dplyr)
library(ggrepel)
library(ggplot2)
library(rafalib)
library(UpSetR)
library(reshape2)
library(ggpubr)
mypar()

## load MPNST and Tumor dataset
load("MPNST.Tumor.RData")
genotypes <- factor(c(rep("Normal", 7), rep("MPNST",7)), levels = c("MPNST", "Normal"))
neural.genes <- c("SOX9", "MSX1", "ZIC1", "ZIC2", "ZIC3",
                  "DLX5", "DLX6", "SNAI2", "ID1", "PAX7", "SOX5", "PAX3")
schwann.genes <- c("MBP", "S100B", "MPZ", "L1CAM", "NCAM2", "GAP43", "SOX2")
neural.genes <- schwann.genes
mat <- mpnst.normal$result$results[mpnst.normal$result$results$gene %in% 
                                     neural.genes,c(7:20,22)]
mat <- mat[order(mat$gene), ]
mat1 <- melt(mat)
mat1$genotypes <- rep(genotypes, each = length(neural.genes))
mat1$value <- log2(mat1$value + 1)
mat1$genotypes <- relevel(mat1$genotypes, ref = "Normal")

## Box plot
p1 <- ggboxplot(mat1, "gene", "value", fill = "genotypes",
                palette = c("#00cccc","#ff4d4d"), notch = TRUE) + theme_grey() +
        facet_grid(.~ gene,  space = "free", scale = "free") + ylab("Log2 Normalized Counts") + 
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
              axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
              axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank(),  
              axis.title.x =  element_blank(), axis.line.x = element_blank(),
              legend.text = element_text(size = 18, face = "bold", color = "black"),
              strip.text = element_text(size = 14, colour = "black", face = "bold"))
p1

## pheatmap
neural.genes <- c("SOX9", "MSX1", "ZIC1", "ZIC2", "ZIC3",
                  "DLX5", "DLX6", "SNAI2", "ID1", "PAX7", "SOX5", "PAX3")
mat2 <- mpnst.normal$result$results[mpnst.normal$result$results$gene %in% 
                                     neural.genes,c(7:21)]
mat2 <- mat2[mat2$Significant == "Up",c(1:14)]
genotypes <- relevel(genotypes, ref = "Normal")
annot1 <- data.frame(genotypes, row.names = colnames(mat2))
pheatmap(mat = log2(mat2+1), cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = F,legend = TRUE, 
         fontsize = 14, scale = "row", fontface="bold", 
         clustering_distance_rows = "correlation", fontsize_row = 18,
         clustering_distance_cols = "euclidean", annotation = annot1,
         filename = "../results/05 Neural Crest Genes/NeuralCrest_MPNST-vs-Normal.png", 
         width = 6, height = 8, border_color = "white")

## Neural Crest for all cell lines
## Clearing Workspace
rm(list = ls())

## libraries used
library(DESeq2)
library(dplyr)
library(genefilter)
library(ggrepel)
library(ggplot2)
library(ggthemes)
library(pheatmap)
library(rafalib)
mypar()

## source all the function
source("combiningFilesinDirectory.R")
source("DESeqCalculation.R")
source("plots2Plot.R")

## files
filePath <- "counts/"
file_list <- list.files(path = filePath, pattern = "(.*)ReadsPerGene.out.tab")
counts.table <- combiningFilesinDirectory(filePath = filePath, 
                                          file_list = file_list, n = 2)
counts.table <- counts.table[-c(1:4),]

## MPNST WT/NF1/Sporadic associated samples for analysis
MPNST007 <- counts.table %>% select(contains("MPNST007")) %>% colnames()
M4970 <- counts.table %>% select(contains("MPNST4970")) %>% colnames()
M181 <- counts.table %>% select(contains("MPNST181")) %>% colnames()
SNF02.2 <- counts.table %>% select(contains("SNF02")) %>% colnames()
MPNST724 <- counts.table %>% select(contains("Sample_MPNST724")) %>% colnames()
STS26T <- counts.table %>% select(contains("Sample_STS26T")) %>% colnames()

## MPNST PRC2 Negative associated samples for analysis
ST88 <- counts.table %>% select(contains("Sample_MPNST-14-JF")) %>% colnames()
S462 <- (counts.table %>% select(contains("Sample_S462-p")) %>% colnames())[4:6]
MPNST642 <- counts.table %>% select(contains("MPNST642")) %>% colnames()
MPNST3813E <- counts.table %>% select(contains("MPNST1003813E")) %>% colnames()
T265 <- counts.table %>% select(contains("MPNSTT265")) %>% colnames()
SNF96.2 <- counts.table %>% select(contains("SNF96")) %>% colnames()

## selecting the PRC2 Negative and WT samples
counts.table <- counts.table[,c(MPNST007, M4970, M181, SNF02.2, MPNST724, STS26T,
                                ST88, S462, MPNST642, MPNST3813E, T265, SNF96.2)]
dim(counts.table)

## DESeq Analysis 
cell.lines <- factor(c(rep("MPNST007",3),rep("MPNST4970",3), rep("MPNST181",3),
                       rep("SNF02.2",3), rep("MPNST724",3), rep("STS26T",3),
                       rep("ST88-14",3), rep("S462",3), rep("MPNST642",3),
                       rep("MPNST3813E",3), rep("T265",3), rep("SNF96.2",3)),
                     levels=c("STS26T","MPNST724","MPNST181",
                              "MPNST007","MPNST4970","SNF02.2",
                              "ST88-14","S462","MPNST642","SNF96.2",
                              "T265","MPNST3813E"))
genotypes <- factor(c(rep("PRC2-WT",18), rep("PRC2-Loss",18)),
                    levels = c("PRC2-WT","PRC2-Loss"))
colData <- data.frame(colnames(counts.table), cell.lines,genotypes)    
dds <- DESeqDataSetFromMatrix(counts.table, colData, design = ~ genotypes)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
idx <- which(rowMeans(dat) > 30)
dat <- dat[idx,]
dds <- dds[idx,]
print(dim(dat))
print(dim(dds))

## neural crest
# neural.genes <- c("SOX9", "MSX1", "ZIC1", "ZIC2", "ZIC3", "DLX5", "DLX6", "SNAI2",
#                   "ID1", "PAX7", "SOX5", "PAX3")
#mat <- dat[rownames(dat) %in% neural.genes,]
cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')
neural.crest <- read.table(file = "neural-crest_genes.txt",
                           header = F, sep = "\t", stringsAsFactors = F)
mat <- dat[rownames(dat) %in% neural.crest$V1, ]
annot1 <- data.frame(genotypes = relevel(colData$genotypes,ref = c("PRC2-WT")),
                     cell.lines = cell.lines, row.names = colnames(mat))
annot_cols <- list(cell.lines = cols, genotypes = c("#4DAF4A", "#E41A1C"))
names(annot_cols$cell.lines) <- levels(cell.lines)
names(annot_cols$genotypes) <- levels(genotypes)
pheatmap(mat = log2(mat+1), cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 15,
         legend = TRUE, scale = "row", fontface="bold", annotation = annot1,
         clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
         clustering_method = "centroid", fontsize = 12, 
         annotation_colors = annot_cols, border_color = "white", 
         filename = "NeuralCrest_PRC2-Loss-vs-WT_new.png", 
         width = 8, height = 12)


