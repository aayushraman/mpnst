## Comparison of between tumor samples
tumor.comp <- function(res.dds, geno, cols, annot, fc, tr.fal, rld){
  message(sum(res.dds$padj < 0.05 & abs(res.dds$log2FoldChange) > log2(fc), na.rm = TRUE))
  results <- as.data.frame(res.dds)
  
  ## Sorted as per adj. P-value
  resSort <- res.dds[order(res.dds$padj),]
  print(head(resSort))
  topGenes <- resSort[1:20,]
  topGenes$genes <- rownames(topGenes)
  
  ## upregulated and downregulated genes
  ind.up = which(results$log2FoldChange > log2(fc) & results$padj < 0.05)
  ind.down = which(results$log2FoldChange < log2(1/fc) & results$padj < 0.05)
  up.reg = results[ind.up, ]
  down.reg = results[ind.down, ]
  results.list = list(results = results, up.reg = up.reg, down.reg = down.reg, counts = dat, rld = rld)
  return(results.list)
}

## MPNST related plots
plots.mpnst <- function(dds, res.dds, cols, fc, lev1, geno.cols, tr.fal, genotypes, rld){
  res.dds$norm.counts <- counts(dds[,cols], normalized=TRUE)
  annotation <- data.frame(Sample.Type = genotypes[cols])
  rownames(annotation) <- colnames(dds[,cols])
  #annotation$Sample.Type <- relevel(annotation$Sample.Type, ref = "Neurofibroma")
  comp <- tumor.comp(res.dds = res.dds, geno = genotypes[cols], cols = cols, 
                     annot = annotation, fc = fc, rld =  rld)
  print(head(comp$results))
  print(head(comp$up.reg))
  print(head(comp$down.reg))
  genes <- rownames(comp$results[comp$results$padj < 0.05 & 
                                   abs(comp$results$log2FoldChange) > log2(fc) & 
                                   !is.na(comp$results$padj),])
  message(sum(comp$results$padj < 0.05 & 
                abs(comp$results$log2FoldChange) > log2(fc), na.rm = TRUE))
  
  ## Annotation Data frame
  df1 <- data.frame(Samples = genotypes[cols])
  geno.cols.assigned <- setNames(geno.cols, unique(levels(df1$Samples))[lev1])
  annot1 <- HeatmapAnnotation(df = df1, col = list(Samples = geno.cols.assigned), show_legend = tr.fal, 
                              annotation_legend_param = list(title = "Genotypes", width = unit(1, "cm"), 
                                                             title_gp = gpar(fontsize = 14,fontface = "bold"), 
                                                             labels_gp = gpar(fontsize = 14, fontface = "bold")))
  mat.scaled <- t(scale(t(assay(comp$rld[genes, cols])),center=TRUE,scale=FALSE))
  plot.all <- Heatmap(mat.scaled, show_row_names = FALSE, show_column_names = TRUE, row_dend_reorder = TRUE, 
                      column_dend_reorder = TRUE, clustering_distance_rows = "pearson",
                      clustering_distance_columns = "euclidean", clustering_method_rows = "complete", 
                      row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                      column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
                      clustering_method_columns = "complete", top_annotation = annot1,
                      heatmap_legend_param = list(title = "Scaled Counts", 
                                                  title_gp = gpar(fontsize = 14, fontface = "bold"), 
                                                  labels_gp = gpar(fontsize = 14, fontface = "bold")))
  ## Top Up and Down Genes
  message("Up = ",sum(comp$up.reg$log2FoldChange > log2(fc) & comp$up.reg$padj < 0.05, na.rm = FALSE))
  genes.up <- rownames(comp$up.reg[order(comp$up.reg$padj),])[1:50]
  plot.up <- Heatmap(mat.scaled[genes.up,], show_row_names = TRUE, show_column_names = TRUE, row_dend_reorder = TRUE, 
                     column_dend_reorder = TRUE, clustering_distance_rows = "pearson",
                     clustering_distance_columns = "euclidean", clustering_method_rows = "complete", 
                     row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                     column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
                     clustering_method_columns = "complete", top_annotation = annot1,
                     heatmap_legend_param = list(title = "Scaled Counts", 
                                                 title_gp = gpar(fontsize = 14, fontface = "bold"), 
                                                 labels_gp = gpar(fontsize = 14, fontface = "bold")))
  message("Down = ", sum(comp$down.reg$log2FoldChange < log2(1/fc) & comp$down.reg$padj < 0.05, na.rm = TRUE))
  genes.down <- rownames(comp$down.reg[order(comp$down.reg$padj),])[1:50]
  plot.down <- Heatmap(mat.scaled[genes.down,], show_row_names = TRUE, show_column_names = TRUE, row_dend_reorder = TRUE, 
                       column_dend_reorder = TRUE, clustering_distance_rows = "pearson",
                       clustering_distance_columns = "euclidean", clustering_method_rows = "complete", 
                       row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                       column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
                       clustering_method_columns = "complete", top_annotation = annot1,
                       heatmap_legend_param = list(title = "Scaled Counts", 
                                                   title_gp = gpar(fontsize = 14, fontface = "bold"), 
                                                   labels_gp = gpar(fontsize = 14, fontface = "bold")))
  
  ## WNT Cat genes
  WNT.Cat <- read.table("../dat-info/WNT-Cat_CommonGenes.txt", header = FALSE, 
                        sep = "\t", stringsAsFactors = FALSE)
  comp$results$Significant <- ifelse(comp$results$log2FoldChange > log2(fc) & 
                                       comp$results$padj < 0.05, "Up",
                                     ifelse(comp$results$log2FoldChange < log2(1/fc) & 
                                              comp$results$padj < 0.05, 
                                            "Down","Not Signif"))
  comp$results$gene <- rownames(comp$results)
  pval <- max(res.dds[which(res.dds$padj < 0.05),"pvalue"])
  dat1 <- subset(comp$results[comp$results$gene %in% WNT.Cat$V1,], 
                -log10(pvalue) > -log10(pval) & abs(log2FoldChange) > log2(fc))
  up.dat1 <- sum(dat1$log2FoldChange > log2(fc), na.rm = TRUE)
  down.dat1 <- sum(dat1$log2FoldChange < log2(fc), na.rm = TRUE)
  cat("Number of DEGs related to WNT/Cadh = ",dim(dat1)[1],"\n")
  cat("Number of Up. WNT/Cadh DEGs = ",up.dat1,"\n")
  cat("Number of Up. WNT/Cadh DEGs = ",down.dat1,"\n")
  plot.Wnt.Cat <- ggplot(comp$results, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) + scale_color_manual(values = c("green", "grey", "red")) + 
    geom_text_repel(data = dat1, aes(label = gene), size = 5, fontface = "bold", color = "black",
                    box.padding = unit(0.35, "lines"), point.padding = unit(0.2, "lines")) + 
    geom_hline(aes(yintercept = -log10(pval)), color="dodgerblue", linetype="dashed") + 
    geom_vline(aes(xintercept = log2(fc)), color="dodgerblue", linetype="dashed") +
    geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue", linetype="dashed") +
    theme_bw(base_size = 16) + theme(axis.title = element_text(size = 22, face = "bold"),
                                     axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                                     axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                                     legend.title = element_text(size = 16, face = "bold", color = "black"),
                                     legend.text = element_text(size = 16, face = "bold", color = "black"),
                                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  plot.Wnt.Cat.withoutlabels <- ggplot(comp$results, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = Significant)) + scale_color_manual(values = c("green", "grey", "red")) + 
      geom_hline(aes(yintercept = -log10(pval)), color="dodgerblue", linetype="dashed") + 
      geom_vline(aes(xintercept = log2(fc)), color="dodgerblue", linetype="dashed") +
      geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue", linetype="dashed") +
      theme_bw(base_size = 16) + theme(axis.title = element_text(size = 22, face = "bold"),
                                       axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                                       axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                                       legend.title = element_text(size = 16, face = "bold", color = "black"),
                                       legend.text = element_text(size = 16, face = "bold", color = "black"),
                                       plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  plot.Wnt.Cat.withoutlabels
  result.list <- list(result = comp, plot.all = plot.all, plot.up = plot.up, 
                      plot.down = plot.down, plot.Wnt.Cat = plot.Wnt.Cat, plot.Wnt.Cat.withoutlabels =  plot.Wnt.Cat.withoutlabels)
  return(result.list)
}

volcano.plot <- function(results, fc, pval){
  ## WNT Cat genes
  WNT.Cat <- read.table("~/Desktop/MPNST/dat-info/WNT-Cat_Patients.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  results$Significant <- ifelse(results$log2FoldChange > log2(fc) & results$padj < 0.05, "Up",
                                ifelse(results$log2FoldChange < log2(1/fc) & results$padj < 0.05, "Down","Not Signif"))
  results$gene <- rownames(results)
  pval <- max(results[which(results$padj < pval),"pvalue"])
  dat1 <- subset(results[results$gene %in% WNT.Cat$V1,], 
                        -log10(pvalue) > -log10(pval) & abs(log2FoldChange) > log2(fc))
  up.dat1 <- sum(dat1$log2FoldChange > log2(fc), na.rm = TRUE)
  down.dat1 <- sum(dat1$log2FoldChange < log2(fc), na.rm = TRUE)
  cat("Number of DEGs related to WNT/Cadh =",dim(dat1)[1],"out of",length(WNT.Cat$V1),"\n")
  cat("Number of Up. WNT/Cadh DEGs = ",up.dat1,"\n")
  cat("Number of Down WNT/Cadh DEGs = ",down.dat1,"\n")
  plot.Wnt.Cat <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) + scale_color_manual(values = c("green", "grey", "red")) + 
    geom_text_repel(data = dat1, aes(label = gene), size = 5, fontface = "bold", color = "black",
                    box.padding = unit(0.35, "lines"), point.padding = unit(0.2, "lines")) + 
    geom_hline(aes(yintercept = -log10(pval)), color="dodgerblue", linetype="dashed") + 
    geom_vline(aes(xintercept = log2(fc)), color="dodgerblue", linetype="dashed") +
    geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue", linetype="dashed") +
    xlab("Log2 Fold-Change") + ylab("-Log10 (P-value)") + theme_bw(base_size = 16) + 
    theme(axis.title = element_text(size = 28, face = "bold"),
          axis.text.x = element_text(size = 26, face = "bold", color = "black"),
          axis.text.y = element_text(size = 26, face = "bold", color = "black"),
          legend.position = "none",
          #legend.title = element_text(size = 20, face = "bold", color = "black"),
          #legend.text = element_text(size = 20, face = "bold", color = "black"),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  plot.Wnt.Cat.withoutlabels <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = Significant)) + scale_color_manual(values = c("green", "grey", "red")) + 
      geom_hline(aes(yintercept = -log10(pval)), color="dodgerblue", linetype="dashed") + 
      geom_vline(aes(xintercept = log2(fc)), color="dodgerblue", linetype="dashed") +
      geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue", linetype="dashed") +
      theme_bw(base_size = 16) + theme(axis.title = element_text(size = 22, face = "bold"),
                                       axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                                       axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                                       legend.title = element_text(size = 16, face = "bold", color = "black"),
                                       legend.text = element_text(size = 16, face = "bold", color = "black"),
                                       plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  return(list(plot.Wnt.Cat = plot.Wnt.Cat, plot.Wnt.Cat.withoutlabels = plot.Wnt.Cat.withoutlabels))
}
