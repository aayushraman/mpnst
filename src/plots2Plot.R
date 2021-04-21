# PCA using ggplots
PCAplot <- function(data, genotypes, conditions, title){
  
  ## Calculating PC components
  pcs = prcomp(t(data), center = TRUE)
  percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2) 
  print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes, shape = genotypes), 
               environment = environment()) + ggtitle(title) + 
          xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + geom_point(size = 8) + 
          theme(legend.text = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 16, colour = "black", face = "bold"),
                plot.title = element_text(size = 0, face ="bold"),
                axis.title = element_text(size = 18, face = "bold"),
                axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}

makeLab = function(x,pc) {
  paste0("PC",pc,": ",x,"% variance")
}

## plotting housekeeping
plotHousekeepingGenes <- function(counts.table, list.Gname, geneNames, genotypes){
  for(i in 1:length(list.Gname)){
    ind = which(rownames(counts.table) == list.Gname[i])  
    plot.dat = data.frame(labels = factor(colnames(counts.table), 
                          levels = colnames(counts.table)),
                          Normalized.Counts = as.vector(as.matrix(counts.table[ind,])),
                          Genotype_Condition = genotypes)
    print(ggplot(plot.dat, aes(x = labels, y = Normalized.Counts, fill = Genotype_Condition)) + 
          geom_bar(stat="identity", width=0.75, position = position_dodge(width=0.5)) + 
          xlab("") + ylab("Log2 Normalized Data") + facet_grid(. ~ Genotype_Condition,  space = "free", scale="free") +
          ggtitle(paste("Barplot for Gene ",geneNames[i],sep="")) + 
            theme(plot.title = element_text(size = 18, face = "bold"),
                  axis.title.y= element_text(size = 18, colour = "black", face = "bold"),
                  axis.text.y = element_text(size = 18, colour = "black", face = "bold"),
                  axis.text.x = element_blank(), legend.position="none", 
                  axis.ticks.x = element_blank(), axis.title.x =  element_blank(),
                  strip.text = element_text(size = 16, colour = "red", face = "bold"))) 
  }
}

## box plot for genes
plotBox_genes <- function(counts.table, list.Gname, geneNames, genotypes){
  for(i in 1:length(list.Gname)){
    ind <- which(rownames(counts.table) == list.Gname[i])  
    plot.dat <- data.frame(labels = factor(colnames(counts.table), levels = colnames(counts.table)),
                          Normalized.Counts = as.vector(as.matrix(counts.table[ind,])), Genotype_Condition = genotypes)
    print(ggplot(plot.dat, aes(x = Genotype_Condition, y = Normalized.Counts, fill = Genotype_Condition)) +  
            geom_boxplot(size = 1) + xlab("") + ylab("Normalized Counts") + geom_point(color = "black", size = 2) + 
            facet_grid(. ~ Genotype_Condition,  space = "free", scale="free") + ylim(c(0,24)) + 
            geom_text_repel(data = plot.dat, aes(label = labels, fontface = "bold", size = 2), segment.color = 'grey50') + 
            ggtitle(paste("Boxplot for Gene ",geneNames[i],sep="")) +
            theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.title.x =  element_blank(),
                  axis.title.y= element_text(size = 18, colour = "black", face = "bold"),
                  axis.text.y = element_text(size = 18, colour = "black", face = "bold"),
                  axis.text.x = element_blank(), legend.position="none", axis.ticks.x = element_blank(), 
                  strip.text = element_text(size = 16, colour = "red", face = "bold"))) 
    cat("Gene: ", list.Gname[i],"\n")
    print(plot.dat)
    cat("\n\n")
  }
}
