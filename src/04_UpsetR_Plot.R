############################################################
# @Author: Ayush T. Raman
# Rai Lab, MD Anderson Cancer Center
# Date: Feb 26th, 2019
#
# Program is used for:
# 1. Common Genes across Celllines and Tumor Datasets
# 2. similar MPNST/src/UpSet_CommonGenes.R 
############################################################

## Clearing Workspace
rm(list = ls())

## libraries used
library(dplyr)
library(ggrepel)
library(ggplot2)
library(rafalib)
library(UpSetR)
mypar()

## loading all the gene expression datasets
load("MPNST.Tumor.RData")
load("results.ST88.RData")
load("results.S462.RData")
load("results.MPNSTE.RData")
load("results.MPNST642.RData")
load("results.T265.RData")
load("results.SNF96.2.RData")

## common genes
genes.cl <- unique(c(results.MPNST642$up.reg$gene, results.MPNST642$down.reg$gene,
                     results.MPNSTE$up.reg$gene, results.MPNSTE$down.reg$gene,
                     results.S462$up.reg$gene, results.S462$down.reg$gene,
                     results.ST88$up.reg$gene, results.ST88$down.reg$gene,
                     results.T265$up.reg$gene, results.T265$down.reg$gene,
                     results.SNF96.2$up.reg$gene, results.SNF96.2$down.reg$gene))
genes.all <- genes.cl[genes.cl %in% c(rownames(mpnst.normal$result$up.reg), 
                                      rownames(mpnst.normal$result$down.reg))]
dat <- data.frame(genes = genes.all, stringsAsFactors = F)

## Common Genes for each type
dat$Tumor <- ifelse(dat$genes %in% rownames(mpnst.normal$result$up.reg), 1, 0)
dat$Tumor.Down <- ifelse(dat$genes %in% rownames(mpnst.normal$result$down.reg), 1, 0)

dat$ST88.14 <- ifelse(dat$genes %in% results.ST88$up.reg$gene, 1, 0)
dat$ST88.14.Down <- ifelse(dat$genes %in% results.ST88$down.reg$gene, 1, 0)

dat$S462 <- ifelse(dat$genes %in% results.S462$up.reg$gene, 1, 0)
dat$S462.Down <- ifelse(dat$genes %in% results.S462$down.reg$gene, 1, 0)

dat$MPNST642 <- ifelse(dat$genes %in% results.MPNST642$up.reg$gene, 1, 0)
dat$MPNST642.Down <- ifelse(dat$genes %in% results.MPNST642$down.reg$gene, 1, 0)

dat$MPNST3813E <- ifelse(dat$genes %in% results.MPNSTE$up.reg$gene, 1, 0)
dat$MPNST3813E.Down <- ifelse(dat$genes %in% results.MPNSTE$down.reg$gene, 1, 0)

dat$T265 <- ifelse(dat$genes %in% results.T265$up.reg$gene, 1, 0)
dat$T265.Down <- ifelse(dat$genes %in% results.T265$down.reg$gene, 1, 0)

dat$SNF96.2 <- ifelse(dat$genes %in% results.SNF96.2$up.reg$gene, 1, 0)
dat$SNF96.2.Down <- ifelse(dat$genes %in% results.SNF96.2$down.reg$gene, 1, 0)

## UpSetR Plot
comp.list1 <- list(list("Tumor", "MPNST642"), list(colnames(dat)[2], "S462"),
                   list("Tumor", "ST88.14"), list("Tumor", "T265"),
                   list("Tumor", "SNF96.2"), list("Tumor", "MPNST3813E"),
                   list("Tumor", "ST88.14", "S462"), list("Tumor", "ST88.14", "MPNST642"),
                   list("Tumor", "S462", "MPNST642"), list("Tumor", "MPNST3813E", "MPNST642"),
                   list("Tumor", "SNF96.2", "MPNST642"),
                   list("Tumor", "MPNST642", "ST88.14", "S462"),
                   list("Tumor", "MPNST642", "ST88.14", "T265"),
                   list("Tumor", "MPNST642", "T265", "S462"),
                   list("Tumor", "MPNST642", "ST88.14", "MPNST3813E"),
                   list("Tumor", "SNF96.2", "S462", "MPNST642"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "S462"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "SNF96.2"),
                   list("Tumor", "MPNST642", "ST88.14", "S462", "MPNST3813E"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "SNF96.2", "MPNST3813E"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "SNF96.2"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "MPNST3813E"),
                   list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "SNF96.2", "MPNST3813E"))
upset(dat[,c(1,2,4,6,8,10,12,14)],nsets = 7, nintersects = 90) 
upset(dat, intersections = comp.list1, main.bar.color = "black", 
      mb.ratio = c(0.5,0.5),text.scale = 2, mainbar.y.label = "Size of Intersection Set", 
      sets.x.label = "Num. of Genes", point.size = 3, line.size = 1.25, 
      group.by = "degree", set_size.angles = 90,
      queries = list(list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "S462"), 
                          color = "dodgerblue1", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265"), 
                          color = "dodgerblue1", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "T265", "S462"), 
                          color = "dodgerblue1", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "MPNST3813E"), 
                          color = "dodgerblue1", active = T),
                     list(query = intersects, params = list("Tumor", "SNF96.2", "S462", "MPNST642"), 
                          color = "dodgerblue1", active = T),

                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "S462"), 
                          color = "mediumorchid3", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "SNF96.2"),
                          color = "mediumorchid3", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "S462", "MPNST3813E"), 
                          color = "mediumorchid3", active = T),
                     
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "SNF96.2", "MPNST3813E"), 
                          color = "mediumseagreen", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "SNF96.2"), 
                          color = "mediumseagreen", active = T),
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "MPNST3813E"), 
                          color = "mediumseagreen", active = T),
                     
                     list(query = intersects, params = list("Tumor", "MPNST642", "ST88.14", "T265", "S462", "SNF96.2", "MPNST3813E"), 
                          color = "red", active = T)))
upset(dat, nsets = 7, nintersects = 90)
