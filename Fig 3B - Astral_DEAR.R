library(IMPRINTS.CETSA.app)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# read RESP data
resp <- openxlsx::read.xlsx("Results/VenetoclaxAstral_RESP_final_summary.xlsx") %>%
  dplyr::filter(kept_for_further_analysis) %>%
  mutate(treatment = sub("Venetoclax", "", treatment)) %>%
  group_by(id, Gene) %>%
  group_modify(~{
    df <- data.frame(combined_pvalue_3h = .x$combined_pvalue[.x$treatment == "3h"],
                     RESP_score_3h = .x$RESP_score[.x$treatment == "3h"],
                     combined_pvalue_12h = .x$combined_pvalue[.x$treatment == "12h"],
                     RESP_score_12h = .x$RESP_score[.x$treatment == "12h"],
                     shift = unique(.x$time_category))
    
    return(df)
  }) %>%
  unique()

# read proteins data
prot3 <- readr::read_tsv("Results/VenetoclaxAstral_proteins_foldchange3h.txt")
prot12 <- readr::read_tsv("Results/VenetoclaxAstral_proteins_foldchange12h.txt")

prot <- full_join(prot3[,-grep("Vehicle|sumUniPeps|sumPSMs|countNum", colnames(prot3))],
                  prot12[,-grep("Vehicle", colnames(prot12))],
                  by = c("id", "description"))
rm(prot3, prot12)

hits <- openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_3h_hits_analysis",
                                       pattern = "summary", full.names = TRUE)) %>%
  rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_hits_analysis",
                                       pattern = "summary", full.names = TRUE))) %>%
  as.data.frame() %>% dplyr::rename(pval = Combinedpval, diff = IS) %>%
  rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_3h_QPHits_analysis",
                                       pattern = "summary", full.names = TRUE)) %>%
          rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_QPHits_analysis",
                                               pattern = "summary", full.names = TRUE))) %>%
          as.data.frame() %>% dplyr::rename(diff = FC, category = criteria_curve) %>% mutate(category = "QP")
  ) %>%
  as.data.frame()

# read DEAR
dear <- openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_3h_DEAR_analysis",
                                       pattern = "summary", full.names = TRUE)) %>%
  rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_DEAR_analysis",
                                       pattern = "summary", full.names = TRUE))) %>%
  as.data.frame() %>% dplyr::select(-criteria_curve) %>%
  dplyr::rename(QPvs37 = FC) %>%
  mutate(treatment = sub(".*-", "", treatment))


dear_hits <- inner_join(dear, unique(hits[,c("id", "Gene", "treatment")]), 
                        by = c("id", "Gene", "treatment"))

dear_hits <- left_join(dear_hits, 
                       resp[,c("id", "Gene", "shift")] %>%
                         dplyr::rename(RESP_shift = shift),
                       by = c("id", "Gene")) %>% 
  mutate(RESP_shift = tidyr::replace_na(RESP_shift, "None"))

### ORA
# getting background
background <- unique(sapply(prot$description, IMPRINTS.CETSA.app:::getGeneName,
                            USE.NAMES = FALSE)) # 7 292 unique proteins

# 3h - GO CC
dearhits_3h_enrich_cc <- clusterProfiler::enrichGO(unique(dear_hits$Gene[which(dear_hits$treatment == "Venetoclax3h")]), 
                                                   ont = "CC", keyType = "SYMBOL",
                                                   OrgDb = "org.Hs.eg.db",  
                                                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                                   minGSSize = 5, maxGSSize = 500,
                                                   universe = background)

dearhits_3h_enrich_cc <- dearhits_3h_enrich_cc@result[which(dearhits_3h_enrich_cc@result$Count >= 5 & 
                                                              dearhits_3h_enrich_cc@result$p.adjust <= 0.05),]
dearhits_3h_enrich_cc <- dearhits_3h_enrich_cc %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  }) %>%
  mutate(treatment = "Venetoclax3h")


# 12h - GO CC
dearhits_12h_enrich_cc <- clusterProfiler::enrichGO(unique(dear_hits$Gene[which(dear_hits$treatment == "Venetoclax12h")]), 
                                                    ont = "CC", keyType = "SYMBOL",
                                                    OrgDb = "org.Hs.eg.db",  
                                                    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                                    minGSSize = 5, maxGSSize = 500,
                                                    universe = background)

dearhits_12h_enrich_cc <- dearhits_12h_enrich_cc@result[which(dearhits_12h_enrich_cc@result$Count >= 5 & 
                                                                dearhits_12h_enrich_cc@result$p.adjust <= 0.05),]
dearhits_12h_enrich_cc <- dearhits_12h_enrich_cc %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  }) %>%
  mutate(treatment = "Venetoclax12h")

# joining both
dearhits_enrich_cc <- as.data.frame(rbind(dearhits_3h_enrich_cc, dearhits_12h_enrich_cc))


### Preparing graph data
relations <- dearhits_enrich_cc[,c("Gene", "treatment", "Description", "p.adjust")] %>%
  dplyr::rename(from = Gene, to = Description, category = treatment, padj = p.adjust) %>% 
  mutate(category = sub("Venetoclax", "DEAR-", category)) %>%
  group_by(category) %>%
  mutate(name = to, padj = -log10(padj)) 

vertex <- dear_hits[, c("Gene", "treatment", "QPvs37", "RESP_shift")] %>%
  mutate(treatment = sub("Venetoclax", "DEAR-", treatment)) %>%
  dplyr::rename(name = Gene, category = treatment) %>%
  mutate(size = 0) %>%
  rbind(unique(dearhits_enrich_cc[,c("Description", "Count", "treatment")]) %>%
          mutate(treatment = sub("Venetoclax", "DEAR-", treatment)) %>%
          dplyr::rename(name = Description, size = Count, category = treatment) %>%
          mutate(category = paste("GO-CC", category), QPvs37 = NA, RESP_shift = NA)) %>%
  mutate(term = ifelse(size != 0, name, NA),
         category = factor(category, levels = c("DEAR-3h", "DEAR-12h", "GO-CC DEAR-3h", "GO-CC DEAR-12h")))

vertex$RESP_shift[which(is.na(vertex$RESP_shift))] <- "None"
vertex$labsize <- 0.4
vertex$labsize[-grep("^DEAR", vertex$category)] <- 2
vertex$size[grep("^DEAR", vertex$category)] <- 2


# plotting graph with relations only at 3 hours
g3 <- igraph::graph_from_data_frame(relations[grep("3h$", relations$category),
                                              c("from", "to", "name", "padj", "category")], 
                                    vertices = vertex[which(grepl("3h$", vertex$category) &
                                                              !is.na(match(vertex$name,
                                                                           unique(c(relations$from[grep("3h$", relations$category)],
                                                                                    relations$to[grep("3h$", relations$category)]))
                                                              )
                                                              )),],
                                    directed = FALSE)

set.seed(123)
g3 <- ggraph(g3, layout = 'kk', weights = padj) + 
  geom_edge_link(aes(width = padj, alpha = padj)) +
  scale_edge_alpha(range = c(0.3,0.75), name = "-log10(p-value)") +
  scale_edge_width(range = c(1,2), name = "-log10(p-value)") +
  geom_node_point(aes(fill = QPvs37, size = size, color = category, 
                      shape = RESP_shift), stroke = 1.8) +
  geom_node_label(aes(x = x, y = y, label = name, fill = QPvs37, size = labsize), 
                  max.time = 2, max.iter = 20000,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 50) +
  scale_size_continuous(breaks = seq(5,30, 10), range = c(3,14),
                        limits = c(0,35), name = "Gene count") +
  scale_shape_manual(values = c("None" = 21, "3h & 12h" = 23, "12h" = 24, "3h" = 25),
                     name = "RESP shift") +
  scale_fill_gradient2(high = "darkorange", low = "blue", name = "QP versus 37°C",
                       na.value = "#cad182", limits = c(-1.5, 1), breaks = c(-1.5,-1,-0.5,0,0.25, 0.5, 1)) +
  scale_color_manual(values = c("DEAR-3h" = "red", "DEAR-12h" = "purple", 
                                "GO-CC DEAR-3h"  = "#8A9039", "GO-CC DEAR-12h" = "#595E0A"), 
                     name = "", na.value = "black",
                     breaks = na.omit(unique(vertex$category))) + 
  guides(color = guide_legend(override.aes = list(size = 5, shape = 21), order = 1),
         fill = guide_legend(override.aes = list(size = 5, shape = 21), order = 2),
         shape = guide_legend(override.aes = list(size = 5), order = 3),
         size = guide_legend(override.aes = list(shape = 21), order = 4)
  ) +
  labs(title = "Gene-concept network from DEAR proteins of the Astral MS Venetoclax dataset at 3 hours",
       subtitle = "adj.p-value cutoff: 0.05 - minimum genes: 5") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 14),
        plot.title = element_text(face = "bold", color = "black", size = 20),
        plot.subtitle = element_text(size = 18))




# plotting graph with relations only at 12 hours
g12 <- igraph::graph_from_data_frame(relations[grep("12h$", relations$category),
                                               c("from", "to", "name", "padj", "category")], 
                                     vertices = vertex[which(grepl("12h$", vertex$category) &
                                                               !is.na(match(vertex$name,
                                                                            unique(c(relations$from[grep("12h$", relations$category)],
                                                                                     relations$to[grep("12h$", relations$category)]))
                                                               )
                                                               )),],
                                     directed = FALSE)

set.seed(123)
g12 <- ggraph(g12, layout = 'kk', weights = padj) + 
  geom_edge_link(aes(width = padj, alpha = padj)) +
  scale_edge_alpha(range = c(0.3,0.75), name = "-log10(p-value)") +
  scale_edge_width(range = c(1,2), name = "-log10(p-value)") +
  geom_node_point(aes(fill = QPvs37, size = size, color = category, 
                      shape = RESP_shift), stroke = 1.8) +
  geom_node_label(aes(x = x, y = y, label = name, fill = QPvs37, size = labsize), 
                  max.time = 2, max.iter = 20000,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  #angle = 0, hjust = 0.5, nudge_y = 0.1, 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 50) +
  scale_size_continuous(breaks = seq(5,30, 10), range = c(3,14),
                        limits = c(0,35), name = "Gene count") +
  scale_shape_manual(values = c("None" = 21, "3h & 12h" = 23, "12h" = 24, "3h" = 25),
                     name = "RESP shift") +
  scale_fill_gradient2(high = "darkorange", low = "blue", name = "QP versus 37°C",
                       na.value = "#cad182", limits = c(-1.5, 1), breaks = c(-1.5,-1,-0.5,0,0.25, 0.5, 1)) +
  scale_color_manual(values = c("DEAR-3h" = "red", "DEAR-12h" = "purple", 
                                "GO-CC DEAR-3h"  = "#8A9039", "GO-CC DEAR-12h" = "#595E0A"), 
                     name = "", na.value = "black",
                     breaks = na.omit(unique(vertex$category))) + 
  guides(color = guide_legend(override.aes = list(size = 5, shape = 21), order = 1),
         fill = guide_legend(override.aes = list(size = 5, shape = 21), order = 2),
         shape = guide_legend(override.aes = list(size = 5), order = 3),
         size = guide_legend(override.aes = list(shape = 21), order = 4)
  ) +
  labs(title = "Gene-concept network from DEAR proteins of the Astral MS Venetoclax dataset at 12 hours",
       subtitle = "adj.p-value cutoff: 0.05 - minimum genes: 5") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 14),
        plot.title = element_text(face = "bold", color = "black", size = 20),
        plot.subtitle = element_text(size = 18))


ggsave("Figures/ORA_GO-CC_network_DEAR_AstralVenetoclax_3h.png",
       plot = g3,  dpi = 300, device = "png", width = 16, height = 10)
ggsave(Figures/"ORA_GO-CC_network_DEAR_AstralVenetoclax_12h.png",
       plot = g12,  dpi = 300, device = "png", width = 16, height = 10)

