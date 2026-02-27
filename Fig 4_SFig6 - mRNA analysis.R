library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)

## 3H
v3h_dmso3h_diff <- openxlsx::read.xlsx("Data/RNAseq/mRNA/V3HvsDMSO3H.xlsx")
v3h_dmso3h_diff$score <- -log10(v3h_dmso3h_diff$pvalue)*sign(v3h_dmso3h_diff$log2FoldChange)
v3h_dmso3h_diff$regulation <- ifelse(v3h_dmso3h_diff$pvalue <= 0.05,
                                     ifelse(v3h_dmso3h_diff$log2FoldChange > 1, "up",
                                            ifelse(v3h_dmso3h_diff$log2FoldChange < -1, "down",
                                                   "no")),
                                     "no")
v3h_dmso3h_diff$regulation[which(is.na(v3h_dmso3h_diff$regulation))] <- "no"

## 12H
v12h_dmso12h_diff <- openxlsx::read.xlsx("Data/RNAseq/mRNA/V12HvsDMSO12H.xlsx")
v12h_dmso12h_diff$score <- -log10(v12h_dmso12h_diff$pvalue)*sign(v12h_dmso12h_diff$log2FoldChange)
v12h_dmso12h_diff$regulation <- ifelse(v12h_dmso12h_diff$pvalue <= 0.05,
                                       ifelse(v12h_dmso12h_diff$log2FoldChange > 1, "up",
                                              ifelse(v12h_dmso12h_diff$log2FoldChange < -1, "down",
                                                     "no")),
                                       "no")
v12h_dmso12h_diff$regulation[which(is.na(v12h_dmso12h_diff$regulation))] <- "no"

## Venn diagramm
mrna_hits_venn <- list("3 hours up" = v3h_dmso3h_diff$gene_id[which(v3h_dmso3h_diff$regulation == "up")],
          "3 hours down" = v3h_dmso3h_diff$gene_id[which(v3h_dmso3h_diff$regulation == "down")],
          "12 hours up" = v12h_dmso12h_diff$gene_id[which(v12h_dmso12h_diff$regulation == "up")],
          "12 hours down" = v12h_dmso12h_diff$gene_id[which(v12h_dmso12h_diff$regulation == "down")])

mrna_hits_venn <- VennDiagram::venn.diagram(mrna_hits_venn, filename = NULL,
                                            disable.logging = TRUE, margin = 0.5,
                                            fill = c("#FF0000", "#008FFF",  "#5400D7", "#00D79F"), lwd = 2, cex = 3,
                                            main = "Kasumi-1 treated with Venetoclax\nRNAseq hits time comparison",
                                            main.fontface = "bold", main.cex = 3, main.pos = c(0.5,0.8),
                                            fontface = "bold", fontfamiliy = "sans",
                                            cat.cex = 2.4, cat.fontface = "bold",
                                            cat.fontfamily = "sans")

png("Figures/RNAseq_Venetoclax_updown_time_VennDiagram.png",
    width=12, height=12, units = "in", res=300)
grid::grid.newpage()
grid::grid.draw(mrna_hits_venn)
dev.off()

### GSEA Tree plot and GSEA curve
d <- godata('org.Hs.eg.db', ont="BP")

### 3H
v3h_dmso3h_gsea <- v3h_dmso3h_diff$score
names(v3h_dmso3h_gsea) <- v3h_dmso3h_diff$gene_id
v3h_dmso3h_gsea <- v3h_dmso3h_gsea[which(!is.na(v3h_dmso3h_gsea))]
v3h_dmso3h_gsea <- sort(v3h_dmso3h_gsea, decreasing = TRUE)
v3h_dmso3h_gsea <- gseGO(v3h_dmso3h_gsea,
                         OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL",
                         pvalueCutoff = 0.01, minGSSize = 10, eps = 1e-16,
                         pAdjustMethod = "BH")

# Treeplot - GSEA pos
x <- v3h_dmso3h_gsea
x@result <- x@result[which(x@result$enrichmentScore > 0),]
x <- pairwise_termsim(x, method = "Lin", semData = d)

t <- treeplot(x, showCategory = 70, label_format = 100,
              cluster.params = list(n = 3, color = c("#F7572E", "#033f63", "#279B5F"), # "#7c9885"), #, "#b5b682", "#B90101"),
                                    label_words_n = 3),
              offset.params = list(tiplab = rel(12), bar_tree = rel(3),
                                   extend = 0.2, hexpand = 0.4))

# setting text appearance of higher classes
t$layers[[3]]$aes_params$size <- 6
t$layers[[3]]$aes_params$fontface <- "bold"
t$layers[[3]]$aes_params$hjust <- -0.02
# setting size of cluster bars
t$layers[[4]]$aes_params$size <- 4
# setting text appearance of GO terms
t$layers[[7]]$aes_params$size <- 4

# add NES score as color
y <- x@result %>%
  dplyr::select(Description, enrichmentScore, NES) %>%
  dplyr::rename(label = Description)
t@data <- left_join(t@data, y, by = "label")

t@layers <- t@layers[-5]
t@layers[[5]]$aes_params$shape <- 21
t@layers[[5]]$aes_params$stroke <- 2
t@layers[[5]]$mapping$fill <- t@data$enrichmentScore

t + scale_fill_gradient2(high = "red", low = "white", name = "Enrichment score",
                         breaks = seq(0,0.8,0.2), limits = c(0,0.8)) +
  scale_color_gradient(name = "p.adjust", trans = "log", high = "#E6A540", low = "#00890B") +
  scale_size_continuous(name = "Gene count", limits = c(0,500),
                        breaks = c(20,60,seq(100,300, 50), 400, 500),
                        range = c(2,16)) +
  labs(title = "GSEA - positively enriched - Venetoclax 3H vs DMSO 3H",
       subtitle = "GO - Biological Process\nBH p-value cutoff: 0.01") +
  theme(plot.title =  element_text(face = "bold", size = 16),
        plot.subtitle =  element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14))

ggsave("Figures/Treeplot_GSEApos_GO-BP_V3HvsDMSO3H.png", device = "png", dpi = 300,
       height = 10, width = 16)


### 12H
v12h_dmso12h_gsea <- v12h_dmso12h_diff$score
names(v12h_dmso12h_gsea) <- v12h_dmso12h_diff$gene_id
v12h_dmso12h_gsea <- v12h_dmso12h_gsea[which(!is.na(v12h_dmso12h_gsea))]
v12h_dmso12h_gsea <- sort(v12h_dmso12h_gsea, decreasing = TRUE)
v12h_dmso12h_gsea <- gseGO(v12h_dmso12h_gsea,
                           OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL",
                           pvalueCutoff = 0.01, minGSSize = 10, eps = 1e-16,
                           pAdjustMethod = "BH")

# Treeplot - GSEA - neg
x <- v12h_dmso12h_gsea
x@result <- x@result[which(x@result$enrichmentScore < 0),]
x <- pairwise_termsim(x, method = "Lin", semData = d)

t <- treeplot(x, showCategory = 50, label_format = 100,
              cluster.params = list(n = 5, color = c("#279B5F",  "#033f63", "#F7572E", "#7c9885", "#B90101"),# "#B90101"),
                                    label_words_n = 4),
              offset.params = list(tiplab = rel(12), bar_tree = rel(3),
                                   extend = 0.2, hexpand = 0.7))

# setting text appearance of higher classes
t$layers[[3]]$aes_params$size <- 6
t$layers[[3]]$aes_params$fontface <- "bold"
t$layers[[3]]$aes_params$hjust <- -0.02
# setting size of cluster bars
t$layers[[4]]$aes_params$size <- 4
# setting text appearance of GO terms
t$layers[[7]]$aes_params$size <- 4

# add NES score as color
y <- x@result %>%
  dplyr::select(Description, enrichmentScore, NES) %>%
  dplyr::rename(label = Description)
t@data <- left_join(t@data, y, by = "label")

t@layers <- t@layers[-5]
t@layers[[5]]$aes_params$shape <- 21
t@layers[[5]]$aes_params$stroke <- 2
t@layers[[5]]$mapping$fill <- t@data$enrichmentScore

t + scale_fill_gradient(high = "white", low = "blue", name = "Enrichment score",
                        breaks = seq(-0.8, 0,0.2), limits = c(0, -0.9)) +
  scale_color_gradient(name = "p.adjust", trans = "log", high = "#E6A540", low = "#00890B") +
  scale_size_continuous(name = "Gene count", limits = c(0,500),
                        breaks = c(20,60,seq(100,300, 50), 400, 500),
                        range = c(2,16)) +
  labs(title = "GSEA - top 50 negatively enriched - Venetoclax 12H vs DMSO 12H",
       subtitle = "GO - Biological Process\nBH p-value cutoff: 0.01") +
  theme(plot.title =  element_text(face = "bold", size = 16),
        plot.subtitle =  element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14))

ggsave("Figures/Treeplot_GSEAneg_GO-BP_V12HvsDMSO12H.png", device = "png", dpi = 300,
       height = 16, width = 16)


# GSEA curves
graph <- enrichplot::gseaplot2(v12h_dmso12h_gsea,
                               geneSetID = which(v12h_dmso12h_gsea@result$enrichmentScore > 0))
graph[[1]]$labels$title <- "GSEA - positively enriched - Venetoclax 12H vs DMSO 12H"
graph[[1]]$labels$subtitle <- "GO - Biological Process\nBH p-value cutoff: 0.01"
graph[[1]] <- graph[[1]] +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.9),
        plot.title =  element_text(face = "bold", size = 16),
        plot.subtitle =  element_text(size = 14),
        legend.text = element_text(size = 14, face = "bold")
  )
graph

ggsave("Figures/GSEAcurve_GSEApos_GO-BP_V12HvsDMSO12H.png", device = "png", dpi = 300,
       height = 6, width = 10)

# data cleaning
rm(d, v3h_dmso3h_gsea, v12h_dmso12h_gsea, x, t, y, graph)

