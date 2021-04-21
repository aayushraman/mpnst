############################################################
# @Author: Ayush T. Raman
# Rai Lab, MD Anderson Cancer Center
# Date: Feb 26th, 2019
#
# Program is used for:
# 1. Box plot for some genes
############################################################

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

## read the file
dat <- read.table(file = "dat-info/box plot gene list.txt", header = T, sep = "\t", 
                  stringsAsFactors = F, quote = "")
dat1 <- dat[,-c(2:7, 44)]

result_folder <- "Box_Plots/"
for(i in 1:nrow(dat1)){
  mat <- melt(dat1[i,], id.vars = "gene")
  mat$value <- log2(mat$value + 1)
  mat$genotypes <- factor(c(rep("PRC2 WT", 18), rep("PRC2 mut", 18)), 
                          levels = c("PRC2 WT", "PRC2 mut"))
  p1 <- ggboxplot(mat, "gene", "value", fill = "genotypes", add = "dotplot",
                  palette = c("#00cccc","#ff4d4d"), notch = TRUE) + theme_grey() +
            facet_grid(.~ gene,  space = "free", scale = "free") + ylab("Log2 Normalized Counts") + 
            theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
                  axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  legend.title = element_blank(), axis.title.x =  element_blank(), 
                  axis.line.x = element_blank(),
                  legend.text = element_text(size = 18, face = "bold", color = "black"),
                  strip.text = element_text(size = 14, colour = "black", face = "bold"))
  plot(p1)
  ggsave(filename = paste(result_folder,"BoxPlot-",dat$gene[i],".tiff",sep = ""),
         plot = p1, width = 10, height = 8, dpi = 600, type = "cairo")
}

mat <- melt(dat1, id.vars = "gene")
mat <- mat[order(mat$gene),]
mat$value <- log2(mat$value + 1)
mat$genotypes <- factor(rep(c(rep("PRC2 WT", 18), rep("PRC2 mut", 18)),12), 
                        levels = c("PRC2 WT", "PRC2 mut"))
p1 <- ggboxplot(mat, "gene", "value", fill = "genotypes", add = "jitter",
                palette = c("#00cccc","#ff4d4d"), notch = TRUE) + theme_grey() +
  facet_grid(.~ gene,  space = "free", scale = "free") + ylab("Log2 Normalized Counts") + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), axis.title.x =  element_blank(), 
        axis.line.x = element_blank(),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        strip.text = element_text(size = 14, colour = "black", face = "bold"))
plot(p1)
ggsave(filename = paste(result_folder,"BoxPlot-All.tiff",sep = ""),
       plot = p1, width = 18, height = 10, dpi = 300, type = "cairo")



