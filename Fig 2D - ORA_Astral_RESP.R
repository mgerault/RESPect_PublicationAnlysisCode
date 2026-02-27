library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

library(IMPRINTS.CETSA.app)
library(ggraph)

### Network summarizing GO term
x <- openxlsx::read.xlsx("VenetoclaxAstral_RESP_final_summary.xlsx") %>%
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
background <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_peptides_postnorm.txt")
# 124 390 peptides
background <- unique(sapply(background$description, IMPRINTS.CETSA.app:::getGeneName,
                            USE.NAMES = FALSE)) # 124 390 peptides from 9 464 gene

# GO BP
x_enrich <- clusterProfiler::enrichGO(x$Gene, ont = "BP", keyType = "SYMBOL",
                                      OrgDb = "org.Hs.eg.db",  
                                      pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                      minGSSize = 5, maxGSSize = 2000,
                                      universe = background)

x_enrich <- x_enrich@result[which(x_enrich@result$Count >= 5 & x_enrich@result$p.adjust <= 0.05),]
# shorter names
x_enrich$Description <- sapply(x_enrich$Description,
                               function(d){
                                 if(nchar(d) > 45){
                                   d_cut <- seq(1, nchar(d), 44)[-1]
                                   if(nchar(d) %in% d_cut)
                                     d_cut <- d_cut[-length(d_cut)]
                                   
                                   d <- strsplit(d, "")[[1]]
                                   d_space <- grep(" ", d)
                                   d_replace <- sapply(d_cut, function(y) d_space[which.min(abs(d_space - y))])
                                   d[d_replace] <- "\n"
                                   d <- paste(d, collapse = "")
                                 };
                                 d
                               }, USE.NAMES = FALSE)


x_enrich <- x_enrich %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  })
colnames(x_enrich)[1:4] <- paste0("GO_", colnames(x_enrich)[1:4])


# GO MF
x_enrich_mf <- clusterProfiler::enrichGO(x$Gene, ont = "MF", keyType = "SYMBOL",
                                         OrgDb = "org.Hs.eg.db",  
                                         pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                         minGSSize = 5, maxGSSize = 2000,
                                         universe = background)

x_enrich_mf <- x_enrich_mf@result[which(x_enrich_mf@result$Count >= 5 & x_enrich_mf@result$p.adjust <= 0.05),]
# no long names
x_enrich_mf <- x_enrich_mf %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  })
colnames(x_enrich_mf)[1:4] <- paste0("GOMF_", colnames(x_enrich_mf)[1:4])


# GO CC
x_enrich_cc <- clusterProfiler::enrichGO(unique(x$Gene), ont = "CC", keyType = "SYMBOL",
                                      OrgDb = "org.Hs.eg.db",  
                                      pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                      minGSSize = 100, maxGSSize = 2000,
                                      universe = background)

x_enrich_cc <- x_enrich_cc@result[which(x_enrich_cc@result$Count >= 5 & x_enrich_cc@result$p.adjust <= 0.05),]
# no long names
x_enrich_cc <- x_enrich_cc %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  })
colnames(x_enrich_cc)[1:4] <- paste0("GOCC_", colnames(x_enrich_cc)[1:4])


# KEGG
k <- clusterProfiler::bitr(unique(sort(x$Gene)),
                           fromType = "SYMBOL", toType = c("ENTREZID"),
                           OrgDb = "org.Hs.eg.db", drop = FALSE)
kbackground <- clusterProfiler::bitr(sort(background),
                                     fromType = "SYMBOL", toType = c("ENTREZID"),
                                     OrgDb = "org.Hs.eg.db", drop = FALSE)
kbackground <- na.omit(kbackground$ENTREZID)

x_enrich_kegg <- clusterProfiler::enrichKEGG(k$ENTREZID, pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                             minGSSize = 5, maxGSSize = 2000,
                                             universe = kbackground)

x_enrich_kegg <- x_enrich_kegg@result[which(x_enrich_kegg@result$Count >= 5 & x_enrich_kegg@result$p.adjust <= 0.05),]
# no long names here
x_enrich_kegg <- x_enrich_kegg %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(ENTREZID = strsplit(.x$geneID, "/")[[1]]) %>%
      left_join(k, by = "ENTREZID") %>%
      dplyr::select(SYMBOL) %>% dplyr::rename(Gene = SYMBOL)
    return(df)
  })
colnames(x_enrich_kegg)[1:4] <- paste0("KEGG_", colnames(x_enrich_kegg)[1:4])

# WikiPathway
x_enrich_wp <- clusterProfiler::enricher(k$ENTREZID,
                                      TERM2GENE = IMPRINTS.CETSA.app:::get_wikipath(FALSE),  
                                      pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                      minGSSize = 5, maxGSSize = 2000,
                                      universe = kbackground)

x_enrich_wp <- x_enrich_wp@result[which(x_enrich_wp@result$Count >= 5 & x_enrich_wp@result$p.adjust <= 0.05),]
x_enrich_wp$Description <- gsub("%.{1,}", "", x_enrich_wp$Description)
# no long names
x_enrich_wp <- x_enrich_wp %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(ENTREZID = strsplit(.x$geneID, "/")[[1]]) %>%
      left_join(k, by = "ENTREZID") %>%
      dplyr::select(SYMBOL) %>% dplyr::rename(Gene = SYMBOL)
    return(df)
  })
colnames(x_enrich_wp)[1:4] <- paste0("WP_", colnames(x_enrich_wp)[1:4])



# Pathway Common (Reactome)
x_enrich_pc <- clusterProfiler::enrichPC(unique(x$Gene),
                                         pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                         minGSSize = 5, maxGSSize = 2000,
                                         universe = background)
x_enrich_pc <- x_enrich_pc@result[which(x_enrich_pc@result$Count >= 5 & x_enrich_pc@result$p.adjust <= 0.05),]
# no long names
x_enrich_pc <- x_enrich_pc %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  group_by(ID, Description, p.adjust, Count) %>%
  group_modify(~{
    df <- data.frame(Gene = strsplit(.x$geneID, "/")[[1]])
    return(df)
  })
colnames(x_enrich_pc)[1:4] <- paste0("PC_", colnames(x_enrich_pc)[1:4])

## save all ORA results
openxlsx::write.xlsx(list("GO-BP" = x_enrich, "GO-MF" = x_enrich_mf, "GO-CC" = x_enrich_cc, 
                          "KEGG" = x_enrich_kegg, "WikiPathway" = x_enrich_wp, "Reactome" = x_enrich_pc),
                     "VenetoclaxAstral_RESP_ORA-Fig2.xlsx")

### Preparing graph data
relations <- x_enrich[,c("Gene", "GO_Description", "GO_p.adjust")] %>%
  dplyr::rename(from = Gene, to = GO_Description, padj = GO_p.adjust) %>% 
  mutate(name = to, padj = -log10(padj)) %>%
  rbind(x_enrich_mf[,c("Gene", "GOMF_Description", "GOMF_p.adjust")] %>%
          dplyr::rename(from = Gene, to = GOMF_Description, padj = GOMF_p.adjust) %>% 
          mutate(name = to, padj = -log10(padj)) %>%
          rbind(x_enrich_cc[,c("Gene", "GOCC_Description", "GOCC_p.adjust")] %>%
                  dplyr::rename(from = Gene, to = GOCC_Description, padj = GOCC_p.adjust) %>% 
                  mutate(name = to, padj = -log10(padj)) %>%
                  rbind(x_enrich_kegg[,c("Gene", "KEGG_Description", "KEGG_p.adjust")] %>%
                          dplyr::rename(from = Gene, to = KEGG_Description, padj = KEGG_p.adjust) %>% 
                          mutate(name = to, padj = -log10(padj)) %>%
                          rbind(x_enrich_wp[,c("Gene", "WP_Description", "WP_p.adjust")] %>%
                                  dplyr::rename(from = Gene, to = WP_Description, padj = WP_p.adjust) %>% 
                                  mutate(name = to, padj = -log10(padj)) %>%
                                  rbind(x_enrich_pc[,c("Gene", "PC_Description", "PC_p.adjust")] %>%
                                          dplyr::rename(from = Gene, to = PC_Description, padj = PC_p.adjust) %>% 
                                          mutate(name = to, padj = -log10(padj)) 
                                        )
                                )
                        )
                )
        )

vertex <- unique(x[,c("Gene", "shift")]) %>%
  dplyr::rename(name = Gene, category = shift) %>%
  mutate(size = 0) %>%
  rbind(unique(x_enrich[,c("GO_Description", "GO_Count")]) %>%
          dplyr::rename(name = GO_Description, size = GO_Count) %>%
          mutate(category = "GO-BP") %>%
          rbind(unique(x_enrich_mf[,c("GOMF_Description", "GOMF_Count")]) %>%
                  dplyr::rename(name = GOMF_Description, size = GOMF_Count) %>%
                  mutate(category = "GO-MF")
                ) %>%
          rbind(unique(x_enrich_cc[,c("GOCC_Description", "GOCC_Count")]) %>%
                  dplyr::rename(name = GOCC_Description, size = GOCC_Count) %>%
                  mutate(category = "GO-CC")
                ) %>%
          rbind(unique(x_enrich_kegg[,c("KEGG_Description", "KEGG_Count")]) %>%
                  dplyr::rename(name = KEGG_Description, size = KEGG_Count) %>%
                  mutate(category = "KEGG")
                ) %>%
          rbind(unique(x_enrich_wp[,c("WP_Description", "WP_Count")]) %>%
                  dplyr::rename(name = WP_Description, size = WP_Count) %>%
                  mutate(category = "WikiPathway")
                ) %>%
          rbind(unique(x_enrich_pc[,c("PC_Description", "PC_Count")]) %>%
                  dplyr::rename(name = PC_Description, size = PC_Count) %>%
                  mutate(category = "Reactome")
                )  
          ) %>% 
  mutate(term = ifelse(size != 0, name, NA),
         category = factor(category, levels = c("3h", "3h & 12h", "12h", 
                                                "GO-BP", "GO-MF", "GO-CC",
                                                "KEGG", "WikiPathway", "Reactome")))

if(any(duplicated(vertex$name))){
  dup_names <- unique(vertex$name[which(duplicated(vertex$name))])
  new_names <- paste(vertex$name[which(!is.na(match(vertex$name, dup_names)))],
                    vertex$category[which(!is.na(match(vertex$name, dup_names)))],
                    sep = " - ")
  
  for(d in dup_names){
    n_names <- vertex$size[vertex$name == d]
    relations$to[which(!is.na(match(relations$to, d)))] <- unlist(mapply(rep, grep(paste0("^", d), new_names, value = TRUE), n_names, 
                                                                         SIMPLIFY = FALSE, USE.NAMES = FALSE))
  }
  
  vertex$term[which(!is.na(match(vertex$name, dup_names)))] <- new_names
  vertex$name[which(!is.na(match(vertex$name, dup_names)))] <- new_names
}
  
#cat_color <- rainbow(length(na.omit(unique(vertex$term))))
#names(cat_color) <- na.omit(unique(vertex$term))
cat_color <- c("RNA splicing, via transesterification\nreactions" = "#FEDC9775",                                    
               "RNA splicing, via transesterification reactions\nwith bulged adenosine as nucleophile" = "#FEDC9775",
               "mRNA splicing, via spliceosome" = "#FEDC9775",  
               "spliceosomal complex" = "#FEDC9775",
               "spliceosomal complex assembly" = "#FEDC9775",
               "Spliceosome" = "#FEDC9775",
               "RNA splicing" = "#FEDC9775",
               "pre-mRNA splicing" = "#FEDC9775",
               
               "RNA processing" = "#28666e75",                                                                      
               "mRNA processing - GO-BP" = "#28666e75",
               "mRNA processing - WikiPathway" = "#28666e75",
               "regulation of mRNA processing" = "#28666e75",
               "mRNA metabolic process" = "#28666e75",    

               "protein-RNA complex assembly" = "#1096A875", 
               "protein-RNA complex organization" = "#1096A875",
               "RNA binding" = "#1096A875",                                                                         
               "mRNA binding" = "#1096A875",                                                                        
               "protein-containing complex binding" = "#1096A875",
               
               
               "regulation of chromosome organization" = "#7c988575",   
               "chromatin organization" = "#7c988575",                                                              
               "chromatin remodeling" = "#7c988575",   
               "chromatin binding" = "#7c988575",
               "chromosome" = "#7c988575",
               "protein localization to organelle" = "#7c988575",                                                   
               "protein localization to chromosome" = "#7c988575",
               "nucleosome organization" = "#7c988575",
               
               "DNA metabolic process" = "#6CA18275",
               "DNA repair" = "#6CA18275",

               "ATP-dependent chromatin remodeling" = "#5E876F75",
               "ATPase complex" = "#5E876F75",
               
               "nucleolus" = "#b5b68275",                                                                           
               "nuclear matrix" = "#b5b68275",                                                                      
               "nuclear body" = "#b5b68275",                                                                        
               "nuclear speck" = "#b5b68275",                                                    
               "nuclear periphery" = "#b5b68275",                                                                   
               "nuclear protein-containing complex" = "#b5b68275",
               
               "rRNA processing" = "#0D72AE75", 
               "rRNA metabolic process" = "#0D72AE75",
               
               "ribonucleoprotein complex biogenesis" = "#033f6375",
               "ribonucleoprotein complex" = "#033f6375", 
               "ribosome biogenesis" = "#033f6375",                                                                 
               "ribosomal large subunit biogenesis" = "#033f6375")

vertex$term[which(is.na(vertex$term))] <- as.character(vertex$category[which(is.na(vertex$term))])

relations <- as.data.frame(rbind(relations,
                                 Reduce(rbind,
                                        sapply(unique(cat_color), 
                                               function(y){
                                                 y <- names(which(cat_color == y))
                                                 
                                                 Reduce(rbind,
                                                        unlist(sapply(y, 
                                                                      function(y1){
                                                                        if(which(y == y1) != length(y))
                                                                          sapply(y[(which(y == y1)+1):length(y)], 
                                                                                 function(y2) data.frame(from = y1, to = y2, padj = NA, name = y1),
                                                                                 simplify = FALSE)
                                                                      }, simplify = FALSE), 
                                                               recursive = FALSE)
                                                 )
                                               }, simplify = FALSE)
                                        )
                                 )
                           )


relations$wgt <- relations$padj
relations$wgt[which(!is.na(relations$wgt))] <- 0.75
relations$wgt[which(is.na(relations$wgt))] <- 0.5

relations <- as.data.frame(rbind(relations,
                                 data.frame(from = c("ribosome biogenesis", "ribosome biogenesis",
                                                     "ATP-dependent chromatin remodeling", "ATPase complex",
                                                     "DNA metabolic process", "DNA repair",
                                                     "RNA processing",
                                                     "mRNA processing - GO-BP", "mRNA processing - WikiPathway", 
                                                     "mRNA metabolic process", "regulation of mRNA processing"),
                                            to = c("rRNA processing", "rRNA metabolic process",
                                                   "chromatin remodeling", "chromatin remodeling",
                                                   "chromatin remodeling", "chromatin remodeling",
                                                   "RNA binding",
                                                   "mRNA binding", "mRNA binding", 
                                                   "mRNA binding", "mRNA binding"),
                                            padj = NA,
                                            name = c("ribosome biogenesis", "ribosome biogenesis",
                                                     "chromatin remodeling", "chromatin remodeling",
                                                     "chromatin remodeling", "chromatin remodeling",
                                                     "RNA processing",
                                                     "mRNA binding", "mRNA binding", 
                                                     "mRNA binding", "mRNA binding"),
                                            wgt = 2
                                            )
                                 )
                           )


# resummarize
relations <- rbind(relations,
                   data.frame(from = "Splicing", 
                              to = unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#FEDC9775"))),
                              padj = NA, name = "Spliceosome", wgt = 0.5)
)

relations <- rbind(relations,
                   data.frame(from = "Nucleus", 
                              to = unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#b5b68275"))),
                              padj = NA, name = "nucleolus", wgt = 0.5)
)

relations <- rbind(relations,
                   data.frame(from = "Chromosome / transcription regulation", 
                              to = unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#7c988575" |
                                                        cat_color[relations$to[which(!is.na(relations$padj))]] == "#6CA18275" |
                                                        cat_color[relations$to[which(!is.na(relations$padj))]] == "#5E876F75"
                              ))),
                              padj = NA, name = "chromosome", wgt = 0.5)
)

relations <- rbind(relations,
                   data.frame(from = "Ribosome biogenesis / processing", 
                              to = unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#0D72AE75" |
                                                        cat_color[relations$to[which(!is.na(relations$padj))]] == "#033f6375" 
                              ))),
                              padj = NA, name = "rRNA processing", wgt = 0.5)
)

relations <- rbind(relations,
                   data.frame(from = "RNA processing / binding", 
                              to = unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#1096A875" |
                                                        cat_color[relations$to[which(!is.na(relations$padj))]] == "#28666e75"
                              ))),
                              padj = NA, name = "RNA binding", wgt = 0.5)
)

relations$to[which(!is.na(relations$padj) &
                         !is.na(match(relations$to, 
                                      unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#FEDC9775")))
                                      ))
                )] <- "Splicing"

relations$to[which(!is.na(relations$padj) &
                     !is.na(match(relations$to, 
                                  unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#b5b68275")))
                     ))
                   )] <- "Nucleus"

relations$to[which(!is.na(relations$padj) &
                     !is.na(match(relations$to, 
                                  unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#7c988575" |
                                                       cat_color[relations$to[which(!is.na(relations$padj))]] == "#6CA18275" |
                                                       cat_color[relations$to[which(!is.na(relations$padj))]] == "#5E876F75"
                                                     )))
                                  ))
                   )] <- "Chromosome / transcription regulation"


relations$to[which(!is.na(relations$padj) &
                     !is.na(match(relations$to, 
                                  unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#0D72AE75" |
                                                       cat_color[relations$to[which(!is.na(relations$padj))]] == "#033f6375" 
                                                     )))
                                  ))
                   )] <- "Ribosome biogenesis / processing"

relations$to[which(!is.na(relations$padj) &
                     !is.na(match(relations$to, 
                                  unique(names(which(cat_color[relations$to[which(!is.na(relations$padj))]] == "#1096A875" |
                                                       cat_color[relations$to[which(!is.na(relations$padj))]] == "#28666e75" 
                                                     )))
                                  ))
                   )] <- "RNA processing / binding"

vertex <- as.data.frame(rbind(vertex,
                              data.frame(name = c("Splicing", "Nucleus", 
                                                  "Chromosome / transcription regulation",
                                                  "Ribosome biogenesis / processing",
                                                  "RNA processing / binding"), 
                                         category = NA,
                                         size = sapply(c("Splicing", "Nucleus", 
                                                         "Chromosome / transcription regulation",
                                                         "Ribosome biogenesis / processing",
                                                         "RNA processing / binding"),
                                                       function(y) length(unique(relations$from[which(!is.na(relations$padj) & relations$to == y)])),
                                                       USE.NAMES = FALSE), 
                                         term = c("Splicing", "Nucleus", 
                                                  "Chromosome / transcription regulation",
                                                  "Ribosome biogenesis / processing",
                                                  "RNA processing / binding")
                                         ) 
                              )
                        )

relations <- relations %>%
  group_by(from, to, wgt) %>% 
  summarise(allterm = paste(name, collapse = "; "),
            padj = paste(padj, collapse = "; "),
            name = name[1])

vertex$labsize <- 0
vertex$labsize[which(is.na(vertex$category))] <- 15
vertex$labsize[grep("[123]", vertex$category)] <- 1
vertex$size[grep("[123]", vertex$category)] <- 5

g <- igraph::graph_from_data_frame(relations, directed = TRUE, 
                                   vertices = vertex[which(!is.na(match(vertex$name,
                                                                        unique(c(relations$from,
                                                                                 relations$to))
                                                                        ))),]
                                   ) 
set.seed(123)
p <- ggraph(g, layout = 'kk', weights = wgt) + 
  geom_edge_link(linewidth = 0.8, alpha = 0.6, aes(color = name), show.legend = FALSE) +
  geom_node_point(aes(color = category, size = size, fill = term),  shape = 21, stroke = 1.8) +
  geom_node_label(aes(x = x, y = y, label = name, fill = term, size = labsize), 
                  max.time = 2, max.iter = 20000,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  #angle = 0, hjust = 0.5, nudge_y = 0.1, 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 50) +
  scale_size_continuous(breaks = seq(5,75, 10), range = c(3,18),
                        limits = c(0,80), name = "Gene count") +
  scale_color_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c", 
                               "GO-BP" = "#2891FF", "GO-MF" = "#28C7FF", "GO-CC" = "#28FFE9",
                               "KEGG" = "#02AA1D", "WikiPathway" = "#5302AA", "Reactome" = "#AA5B02"), 
                    name = "Category", na.value = "black",
                    breaks = na.omit(unique(vertex$category))) +
  scale_fill_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c", 
                               "Splicing" = "#fedc97", "Nucleus" = "#b5b682", 
                               "Chromosome / transcription regulation" = "#7c9885",
                               "Ribosome biogenesis / processing" = "#033f63", 
                               "RNA processing / binding" = "#28666e",
                               cat_color), 
                     name = "GO term", na.value = "black",
                     breaks = na.omit(unique(vertex$term))) + 
  scale_edge_color_manual(values = cat_color) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         fill = "none") +
  labs(title = "Gene-concept network from RESP proteins of the Astral MS Venetoclax dataset",
       subtitle = "adj.p-value cutoff: 0.05 - minimum genes: 5") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14),
        plot.title = element_text(face = "bold", color = "black", size = 20),
        plot.subtitle = element_text(size = 18))


# adding proteins without any assigned terms
knitr::kable(x[!(unique(x$Gene) %in% relations$from),c("Gene", "shift")])
# |Gene     |shift    |
# |:--------|:--------|
# |HSPA4L   |3h & 12h |
# |ARHGAP25 |3h & 12h |
# |PRDX4    |3h & 12h |
# |GIT2     |3h & 12h |
# |GAPVD1   |3h & 12h |
# |KIAA1143 |12h      |
# |SH3KBP1  |3h & 12h |
# |NRBF2    |12h      |

p <- p +
  geom_point(data = data.frame(Gene = c("HSPA4L", "ARHGAP25", "PRDX4", "GIT2", "GAPVD1", "SH3KBP1"),
                               x = c(-6.5, -6.1, 5.5, 5, 5.7, -5.8), y = c(1.2, 0.9, 5.5, 4.9, 2, 4)),
             aes(x, y), shape = 21, stroke = 1.8, size = 6.5,
             color = "#e7487c", fill = "#e7487c") +
  ggrepel::geom_label_repel(data = data.frame(Gene = c("HSPA4L", "ARHGAP25", "PRDX4", "GIT2", "GAPVD1", "SH3KBP1"), 
                                              x = c(-6.5, -6.1, 5.5, 5, 5.7, -5.8), y = c(1.2, 0.9, 5.5, 4.9, 2, 4)),
                            aes(x, y, label = Gene), size = 4.75,
                            color = "black", fill = "#e7487c", fontface = "bold",
                            min.segment.length = 0, max.overlaps = 50)

p <- p +
  geom_point(data = data.frame(Gene = c("NRBF2", "KIAA1143"), 
                               x = c(-0.2,0.2), y = c(-5.7,-5)),
             aes(x, y), shape = 21, stroke = 1.8, size = 6.5,
             color = "#d090f8", fill = "#d090f8") +
  ggrepel::geom_label_repel(data = data.frame(Gene = c("NRBF2", "KIAA1143"), 
                                              x = c(-0.2,0.2), y = c(-5.7,-5)),
                            aes(x, y, label = Gene), size = 4.75,
                            color = "black", fill = "#d090f8", fontface = "bold",
                            min.segment.length = 0, max.overlaps = 50)


ggsave("Figures/Fig2D-ORA_summary_network_RESP_AstralVenetoclax.png",
       plot = p,  dpi = 300, device = "png", width = 26, height = 16)


