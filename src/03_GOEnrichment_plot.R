rm(list = ls())
library(ggplot2)

## Enrichment Analysis
panther.pathways <- read.table(file = "PathwayEnrichment_analysis.txt", sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
colnames(panther.pathways) <- c("pathways", "number", "overUnder", "FDR")
panther.pathways$`Gene Ratio in %` <- panther.pathways$number/14728 * 100
panther.pathways$Regulation <- ifelse(test = panther.pathways$overUnder == "+", "Up.reg", "Down.reg")
panther.pathways$Regulation <- factor(panther.pathways$Regulation, levels = c("Up.reg", "Down.reg"))
panther.pathways$pathways <- gsub(pattern = " \\(P\\d+\\)", replacement = "", panther.pathways$pathways)

## plot
ggplot(data = panther.pathways, aes(y = reorder(pathways, -(FDR)), x = -10*log10(FDR))) + theme_bw(base_size = 16) +
  geom_point(aes(color = Regulation, size = `Gene Ratio in %`)) + xlab("-10*log10(FDR)") + ylab("Panther Pathways") +
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        legend.title = element_text(size = 14, face = "bold", color = "black"), 
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        axis.title = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
