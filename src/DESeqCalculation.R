## DESeq Run
DESeqCalculation <- function(seqcData, genotype, dataset, genes, gene.names, contrast){
  dds <- DESeqDataSetFromMatrix(seqcData, genotype, design = ~ genotypes)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized=TRUE)
  idx <- which(rowMeans(dat) > 30) # 10 for EZH2 KO and MPNST Cell line otherwise 30
  dat <- dat[idx,]
  print(dim(dat))
  dds <- dds[idx,]
  
  ## size factor vs read counts
  sizeFactorRatio <- sizeFactors(dds)
  readCounts <- colSums(counts(dds))
  p1 <- qplot(x = sizeFactors(dds), y = colSums(counts(dds))/1e6) + 
        geom_point(aes(colour = dds$genotypes), size = 5) + theme_classic() +
        geom_smooth(method = "lm", se = TRUE, colour = "grey30") + 
        geom_text_repel(aes(label = colnames(assay(dds))), size = 5, 
                        fontface = "bold", color = "black", 
                        box.padding = unit(0.35, "lines"), point.padding = unit(0.2, "lines")) + 
        xlab("Size Factor") + ylab("Number of Aligned Reads (in million)") +
        theme(axis.title = element_text(size = 16, face = "bold"),
              axis.text.x = element_text(size = 16, face = "bold", color = "black"),
              axis.text.y = element_text(size = 16, face = "bold", color = "black"),
              plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), legend.position = c(.9, .1),
              legend.title = element_blank(), 
              legend.text = element_text(size = 14, face = "bold", color = "black"))
  print(p1)

  ## variability between the samples
  par(mfrow=c(2,3))
  for(i in 1:3){
    for(j in 2:3){
      if((i != j) & (j > i)){
        plot(log2(counts(dds, normalized=TRUE)[,c(i,j)] + 1), pch=16, cex=0.3, font.lab=2)
        abline(0,1, col="red")
      }
    }
  }
  for(i in 4:6){
    for(j in 5:6){
      if((i != j) & (j > i)){
        plot(log2(counts(dds, normalized=TRUE)[,c(i,j)] + 1), pch=16, cex=0.3, font.lab=2)
        abline(0,1, col="red")
      }
    }
  }

  ## Concentration of House keeping genes
  plotHousekeepingGenes(counts.table = log2(dat), list.Gname = genes,
                        genotypes = genotypes, geneNames = gene.names)
  plotBox_genes(counts.table = log2(dat), list.Gname = genes,
                genotypes = genotypes, geneNames = gene.names)
  
  ## PCA Plot
  rld.dds <- assay(rlog(dds, blind=FALSE))
  print(head(rld.dds))
  p2 <- PCAplot(data = rld.dds, genotypes = colData$genotypes, 
                conditions = colData$genotypes, title = "")
  p2 <- p2 + theme_classic()  
  p2 <- p2 + geom_text_repel(aes(label = colData$samples), size = 6, fontface = "bold", 
                                  color = "black", box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.2, "lines")) +
          theme(axis.title = element_text(size = 22, face = "bold"), 
                legend.title = element_blank(),
                axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                legend.text = element_text(size = 16, face = "bold", color = "black"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  print(p2)

  ## Expression Test
  dds <- DESeq(dds, betaPrior = TRUE)
  res.dds <- results(dds, contrast = contrast)
  print(res.dds)
  res.dds$norm.counts <- counts(dds, normalized=TRUE)
  message(sum(res.dds$padj < 0.05 & abs(res.dds$log2FoldChange) > log2(1.5), na.rm = TRUE))
  results <- as.data.frame(res.dds)
  
  ## Sorted as per adj. P-value
  resSort <- res.dds[order(res.dds$padj),]
  print(head(resSort))
  topGenes <- resSort[1:20,]
  topGenes$genes <- rownames(topGenes)
  
  ## Histogram and MA Plot with top 10 genes
  par(mfrow=c(1,1))
  hist(res.dds$pvalue[res.dds$baseMean > 1], breaks=0:20/20, col="grey50", 
       border="white", main="Histogram of p-values with baseMean > 1")
  plotMA(resSort, main = "MA Plot")
  for(i in 1:nrow(topGenes)){
    with(topGenes[i, ],{
      points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
      text(baseMean, log2FoldChange, genes, pos=2, col="dodgerblue")})
  }
  
  ## Volcano Plot
  results <- results %>% tibble::rownames_to_column(var = "gene")
  results <- results[which(!is.na(results$padj)),]
  results$Significant <- ifelse(results$log2FoldChange > log2(1.5) & results$padj < 0.05, "Up",
                                ifelse(results$log2FoldChange < log2(1/1.5) & 
                                         results$padj < 0.05, "Down","Not Signif"))
  pval.sig <- max(results[which(results$padj < 0.05),"pvalue"])
  p3 <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) + 
        geom_point(aes(color = Significant)) + 
        scale_color_manual(values = c("green", "grey", "red")) + theme_classic() +
        xlab("Log2 Fold Change") + ylab("-Log10 P-value") +
        geom_hline(aes(yintercept = -log10(pval.sig)), color="dodgerblue", linetype="dashed") + 
        geom_vline(aes(xintercept = log2(1.5)), color="dodgerblue", linetype="dashed") +
        geom_vline(aes(xintercept = log2(1/1.5)), color="dodgerblue", linetype="dashed") +
        # geom_text_repel(data = subset(results, padj < pval.sig & abs(log2FoldChange) > 1.5),
        #                 aes(label = gene), size = 5, fontface = "bold", color = "black", 
        #                 box.padding = unit(0.35, "lines"), point.padding = unit(0.2, "lines")) +
        theme(axis.title = element_text(size = 22, face = "bold"),
              axis.text.x = element_text(size = 22, face = "bold", color = "black"),
              axis.text.y = element_text(size = 22, face = "bold", color = "black"),
              legend.title = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
              legend.text = element_text(size = 16, face = "bold", color = "black"))
  print(p3)
  
  ## MDS plot
  mdsDist <- cmdscale(d = dist(t(dat)), eig = TRUE, k = 2)
  mdsDist <- cbind(genotype, x = mdsDist$points[,1], y = mdsDist$points[,2])
  p4 <- ggplot(mdsDist, aes(x = x, y = y, color = genotypes)) + theme_classic() +
          xlab("MDS Coordinate 1") + ylab("MDS Coordinate 2") + geom_point(size = 8) + 
          geom_text_repel(aes(label = rownames(mdsDist)), fontface = 'bold', color = 'black',
                          box.padding = unit(0.25, "lines"),
                          point.padding = unit(0.5, "lines")) + 
          theme(legend.text = element_text(size = 14, face = "bold"),
                legend.title = element_text(size = 14, colour = "black", face = "bold"),
                plot.title = element_blank(), axis.title = element_text(size = 18, face = "bold"),
                axis.text.x = element_text(size = 18, face = "bold", color = "black"),
                axis.text.y = element_text(size = 18, face = "bold", color = "black"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  print(p4)
  
  ## legend
  # names(mat_colors$Cell.Lines) <- levels(Genotype)
  # names(mat_colors$Genotype) <- levels(Genotype)

  ## normalized values for the heatmaps
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  mat <- rld.dds[topVarGenes[1:50], ]
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  annot1 <- data.frame(Genotype = relevel(colData$genotypes,ref = c("NSC")),  
                       row.names = colnames(mat))
  p2.top <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 12,
                     legend = TRUE, filename = NA, scale = "row", fontface="bold", 
                     clustering_distance_rows = "correlation", 
                     clustering_distance_cols = "euclidean", annotation = annot1)
  
  ## upregulated and downregulated genes
  ind.up <- which(results$log2FoldChange > log2(1.5) & results$padj < 0.05)
  ind.down <- which(results$log2FoldChange < log2(1/1.5) & results$padj < 0.05)
  up.reg <- results[ind.up, ]
  down.reg <- results[ind.down, ]
  
  ## upregulated heatmaps
  up.rld.dds <- rld.dds[rownames(rld.dds) %in% up.reg$gene,]
  up.topVarGenes <- order(rowVars(up.rld.dds),decreasing=TRUE)
  up.mat <- up.rld.dds[up.topVarGenes[1:50], ]
  p2.up <- pheatmap(up.mat, cluster_rows = TRUE, cluster_cols = TRUE, 
                    show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 12, 
                    legend = FALSE, filename = NA, fontsize_col = 10, scale = "row", 
                    fontface="bold", clustering_distance_rows = "correlation", 
                    clustering_distance_cols = "euclidean", annotation = annot1)

  ## downregulated heatmaps
  down.rld.dds <- rld.dds[rownames(rld.dds) %in% down.reg$gene,]
  down.topVarGenes <- order(rowVars(down.rld.dds),decreasing=TRUE)
  down.mat <- down.rld.dds[down.topVarGenes[1:50], ]
  p3.down <- pheatmap(down.mat, cluster_rows = TRUE, cluster_cols = TRUE, 
                      show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 12, 
                      legend = FALSE, filename = NA, scale = "row", fontface="bold",
                      clustering_distance_rows = "correlation", 
                      clustering_distance_cols = "euclidean", annotation = annot1)
  results.list <- list(results = results, up.reg = up.reg, down.reg = down.reg, counts = dat,
                      plots.p1 = p1, plots.p2 = p2, plots.p3 = p3, plots.p4 = p4)
  return(results.list)
}