###############################################################################################
# @Author: Ayush T. Raman
# Rai Lab, Genome Center, MD Anderson
# Date: August 10th, 2018
#
# Program is used for:
# 1. Clustering Graphs for Histone Marks between various groups
###############################################################################################

rm(list =ls())
setwd("/Volumes/My_Passport/ChromExploreR_all_MPNST_cell_lines/src")

## Library
library(CancerSubtypes)
library(tidyverse)
library(gtools)
library(ggplot2)
library(pheatmap)

## user inputs
sample <- "../dat/sampleNames_MPNST-celllines.txt"
col.sample <- c("Samples", "PRC2 Status")

## sample type labels for plot
colData <- read.table(file = sample, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
colnames(colData) <- col.sample
#chromHMM_folder <- "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E"
bps <- "1000bps"
result_folder <- "../results_new/"

#state_num <- 1

## combined matrix state file from chromStatesMatrixforClustering_particularStates.pl
comp.grp <- function(state_num, colData, file){
  dat <- read.table(file, header = TRUE,
                    sep = "\t",quote = "", row.names = 1, na.strings = FALSE, stringsAsFactors = FALSE)
  idx <- mixedsort(rownames(dat))
  dat <- dat[idx, ]
  print(dim(dat))
  dat.state1 <- FSbyMAD(dat, cut.type="topk",value=10000)
  print(dim(dat.state1))
  genotypes <- factor(colData$`PRC2 Status`, levels = c("NSC","PRC2_Loss","PRC2_WT"))
  sample.name <- colData$Samples
  
  ## binarize the dataset and MDS Plot
  dat.state2 <- ifelse(dat.state1 > 0, 1, 0)
  mds.plot <- cmdscale(d = dist(t(dat.state2), method = "binary"), k = 2)
  print(ggplot(data = data.frame(mds.plot), aes(mds.plot[,1],mds.plot[,2], color = genotypes)) +
          xlab("Comp 1") + ylab("Comp 2") + geom_point(size = 4) + 
          ggrepel::geom_text_repel(aes(label = sample.name)) + theme_classic() +
          theme(legend.text = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 16, colour = "black", face = "bold"),
                plot.title = element_text(size = 0, face ="bold"),
                axis.title = element_text(size = 18, face = "bold"),
                axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  
  ## Wilcoxon Rank Sum Test
  idx1 <- which(colData$`PRC2 Status` == "PRC2_Loss")
  idx2 <- which(colData$`PRC2 Status` != "PRC2_Loss")
  dat.state3 <- data.frame(dat.state1)
  dat.state3$p.values <- as.vector(apply(as.matrix(dat.state1), 1, function(r){
    wilcox.test(r[idx1], r[idx2], exact = FALSE, correct = FALSE)$p.value}))
  hist(as.vector(dat.state3$p.values))
  #dat.state3 <- dat.state3 %>% rownames_to_column(var = "chr_regions")
  dat.state3 <- dat.state3[dat.state3$p.values < 0.05,]
  cat("Total Number of DERs =",dim(dat.state3)[1])
  
  ## heatmap -- 2
  tot.cols <- dim(colData)[1]
  annot1 <- data.frame(genotypes, row.names = colnames(dat.state1)[1:tot.cols])
  annot1$genotypes <- factor(annot1$genotypes, levels = c("NSC","PRC2_Loss","PRC2_WT"))
  pheatmap(dat.state1, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
           show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA,
           fontsize_col = 10, scale = "row", fontface="bold", clustering_distance_rows = "correlation", 
           clustering_distance_cols = "euclidean", annotation = annot1)
  
  ## heatmap -- 2
  #annot1 <- data.frame(genotypes, row.names = colnames(dat.state3)[1:tot.cols])
  #annot1$genotypes <- factor(annot1$genotypes, levels = c("NSC","PRC2_Loss","PRC2_WT"))
  pheatmap(dat.state3[,c(1:tot.cols)], cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
           show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA,
           fontsize_col = 10, scale = "row", fontface="bold", clustering_distance_rows = "correlation", 
           clustering_distance_cols = "euclidean", annotation = annot1)
  
  ## files for DER output and GREAT output
  file1 <- paste(result_folder,"E",state_num,"_DERs.txt",sep = "")
  file2 <- paste(result_folder,"E",state_num,"_GREAT.txt",sep = "")
  bedfile <- data.frame()
  for(j in 1:nrow(dat.state3)){
    bedfile[j, 1] <- noquote(strsplit(rownames(dat.state3)[j], split = ":|-")[[1]])[1]
    bedfile[j, 2] <- noquote(strsplit(rownames(dat.state3)[j], split = ":|-")[[1]])[2]
    bedfile[j, 3] <- noquote(strsplit(rownames(dat.state3)[j], split = ":|-")[[1]])[3]
    bedfile[j, 4] <- paste("ChromState",state_num,"-",j,sep = "")
    bedfile[j, 5] <- sum(dat.state3[j,c(1:tot.cols)])
    bedfile[j, 6] <- sum(dat.state3[j,idx1])
    bedfile[j, 7] <- sum(dat.state3[j,idx2])
    bedfile[j, 8] <- mean(as.numeric(dat.state3[j,idx1]))
    bedfile[j, 9] <- mean(as.numeric(dat.state3[j,idx2]))
  }
  colnames(bedfile) = c("chrom", "chromStart", "chromEnd", "name", "Total Frequency", "# Samples (PRC2-Loss)", "# Samples (PRC2-WT)",
                        "Mean Occurence (PRC2-Loss)", "Mean Occurence (PRC2-Loss)")
  write.table(x = bedfile, file = file1, sep="\t", quote=F, row.names = FALSE, col.names = TRUE)
  write.table(x = bedfile[,c(1:3)], file = file2, sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
}

bed1 <- comp.grp(state_num = 1, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E1-1000bps.txt")
bed2 <- comp.grp(state_num = 2, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E2-1000bps.txt")
bed3 <- comp.grp(state_num = 3, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E3-1000bps.txt")
bed4 <- comp.grp(state_num = 4, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E4-1000bps.txt")
bed5 <- comp.grp(state_num = 5, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E5-1000bps.txt")
bed6 <- comp.grp(state_num = 6, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E6-1000bps.txt")
bed7 <- comp.grp(state_num = 7, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E7-1000bps.txt")
bed8 <- comp.grp(state_num = 8, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E8-1000bps.txt")
bed9 <- comp.grp(state_num = 9, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E9-1000bps.txt")
bed10 <- comp.grp(state_num = 10, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E10-1000bps.txt")
bed11 <- comp.grp(state_num = 11, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E11-1000bps.txt")
bed12 <- comp.grp(state_num = 12, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E12-1000bps.txt")
bed13 <- comp.grp(state_num = 12, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E13-1000bps.txt")
bed14 <- comp.grp(state_num = 12, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E14-1000bps.txt")
bed15 <- comp.grp(state_num = 12, colData = colData, file = "/Volumes/My_Passport/ChromHMM_All_MPNST_cell_lines/ChromHMM_all_lines/15_state/CombinedMatrix-E15-1000bps.txt")

