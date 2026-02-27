library(IMPRINTS.CETSA.app)
library(STRINGdb)
library(ggraph)
library(ggnewscale)

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
                                       pattern = "summary", full.names = TRUE)
                            ) %>%
  rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_hits_analysis",
                                       pattern = "summary", full.names = TRUE)
                            )
        ) %>%
  as.data.frame() %>% dplyr::rename(pval = Combinedpval, diff = IS) %>%
  rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_3h_QPHits_analysis",
                                       pattern = "summary", full.names = TRUE)
                            ) %>%
          rbind(openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_QPHits_analysis",
                                               pattern = "summary", full.names = TRUE)
                                    )
                ) %>%
          as.data.frame() %>% dplyr::rename(diff = FC, category = criteria_curve) %>% mutate(category = "QP")
        ) %>%
  as.data.frame()


# get interactors from each resp using STRINGdb --> top N neighbours
if(!file.exists("STRING_data")){
  dir.create("STRING_data")
}
string_db <- STRINGdb::STRINGdb$new(species=9606,               #ID 9606 correspond to human
                                    score_threshold=200,
                                    input_directory=  file.path(getwd(), "STRING_data"),
                                    version = ifelse(packageVersion("STRINGdb") >= '2.12.0', "12.0", "11.5"))
string_db_proteins <- string_db$get_proteins()

string_idbase <- string_db$map(resp %>% filter(shift == "3h & 12h") %>% 
                                 mutate(id2 = sub(";.*", "", id)) %>%
                                 as.data.frame(), "id2", removeUnmappedRows = FALSE)
string_idbase$id[which(duplicated(string_idbase$id))] # P43243 is duplicated, aka MATR3 --> took isoform MATR3-2
string_db_proteins[which(!is.na(match(string_db_proteins$protein_external_id, string_idbase$STRING_id[string_idbase$id == "P43243"]))),]
# 9606.ENSP00000354346 --> MATR3-2 --> remove it:
string_idbase <- string_idbase[-which(string_idbase$STRING_id == "9606.ENSP00000354346"),]

nmin_interactors <- 50
string_id <- cbind(string_idbase,
                   as.data.frame(Reduce(rbind,
                                        sapply(string_idbase$STRING_id,
                                               function(x){
                                                 gn <- string_idbase$Gene[which(string_idbase$STRING_id == x)] 
                                                 
                                                 neigh <- string_db$get_neighbors(x)
                                                 interactors <- string_db$get_interactions(c(x, neigh))
                                                 interactors <- unique(interactors[which(interactors$from == x | interactors$to == x),])
                                                 interactors <- interactors[na.omit(order(interactors$combined_score, decreasing = TRUE)[1:nmin_interactors]),]
                                                 
                                                 interactors$from_gene <- sapply(interactors$from,
                                                                                 function(g) string_db_proteins$preferred_name[which(string_db_proteins$protein_external_id == g)],
                                                                                 USE.NAMES = FALSE)
                                                 interactors$to_gene <- sapply(interactors$to,
                                                                               function(g) string_db_proteins$preferred_name[which(string_db_proteins$protein_external_id == g)],
                                                                               USE.NAMES = FALSE)
                                                 
                                                 data.frame(top_interactors_string_id = paste(setdiff(unique(c(interactors$from, interactors$to)), x), collapse = "/"),
                                                            top_interactors_gene = paste(setdiff(unique(c(interactors$from_gene, interactors$to_gene)), gn), collapse = "/")
                                                            )
                                               }, simplify = FALSE, USE.NAMES = FALSE)
                                        )
                                 )
                   )

# map them to our cetsa data
string_id$top_interactors_id_cetsa <- sapply(string_id$top_interactors_gene,
                                             function(x){
                                               x <- strsplit(x, "/")[[1]]
                                               #x <- prot$id[unlist(sapply(x, 
                                               #                          function(g) grep(paste0("GN=", g, " "), prot$description))
                                               #                    )];
                                               x <- hits$id[which(!is.na(match(hits$Gene, x)))]
                                               paste(unique(x), collapse = "/")
                                             }, USE.NAMES = FALSE)

string_id$top_interactors_nid_cetsa <- sapply(strsplit(string_id$top_interactors_id_cetsa, "/"), length)

### get all mapped proteins and get their STRING NETWORK
resp_string <- unique(c(unlist(strsplit(string_id$top_interactors_id_cetsa, "/")), string_id$id))
resp_string <- IMPRINTS.CETSA::imprints_average(prot[which(!is.na(match(prot$id, resp_string))),],
                                                savefile = FALSE) 
resp_string$Gene <- sapply(resp_string$description, IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE)
resp_string$description <- sapply(resp_string$description, IMPRINTS.CETSA.app:::getProteinName, USE.NAMES = FALSE)
colnames(resp_string) <- sub("36C", "QP", colnames(resp_string))

# compute I-score (as in imprints_IS)
resp_string$IS_Ventoclax3h <- apply(resp_string[,grep("^\\d{2}C_Venetoclax3h", colnames(resp_string))], 1,
      function(d){
        d = as.numeric(d)
        stabilization = sign(mean(d, na.rm = TRUE))
        d = abs(d)
        
        crit = FALSE
        if(sum(is.na(d)) == length(d)){
          crit = TRUE
        }
        
        if(crit){
          res = NA
        }
        else{
          ord = order(d, decreasing = TRUE)
          d = d[ord]
          dat = data.frame(x = 1:3,
                           y = d)
          wght = rep(0.1, length(d))
          wght[1:2] <- 0.5 # setting greater weight to top 2 FC
          res = lm(y ~ x, data = dat, weights = wght)$coefficients[1]
        };
        res = abs(res)*stabilization
      })

resp_string$IS_Ventoclax12h <- apply(resp_string[,grep("^\\d{2}C_Venetoclax12h", colnames(resp_string))], 1,
                                    function(d){
                                      d = as.numeric(d)
                                      stabilization = sign(mean(d, na.rm = TRUE))
                                      d = abs(d)
                                      
                                      crit = FALSE
                                      if(sum(is.na(d)) == length(d)){
                                        crit = TRUE
                                      }
                                      
                                      if(crit){
                                        res = NA
                                      }
                                      else{
                                        ord = order(d, decreasing = TRUE)
                                        d = d[ord]
                                        dat = data.frame(x = 1:3,
                                                         y = d)
                                        wght = rep(0.1, length(d))
                                        wght[1:2] <- 0.5 # setting greater weight to top 2 FC
                                        res = lm(y ~ x, data = dat, weights = wght)$coefficients[1]
                                      };
                                      res = abs(res)*stabilization
                                    })

resp_string <- resp_string[,-grep("^\\d{2}C_", colnames(resp_string))]

# string mapping
resp_string <- string_db$map(resp_string, "Gene", removeUnmappedRows = FALSE)

resp_string_interactions <- string_db$get_interactions(resp_string$STRING_id)
resp_string_interactions <- unique(resp_string_interactions[resp_string_interactions$combined_score >= 400,])

resp_string2 <- resp_string
resp_string_interactions2 <- resp_string_interactions

resp_string_interactions2$from <- sapply(resp_string_interactions2$from,
                                        function(x) resp_string2$Gene[resp_string2$STRING_id == x],
                                        USE.NAMES = FALSE)
resp_string_interactions2$to <- sapply(resp_string_interactions2$to,
                                        function(x) resp_string2$Gene[resp_string2$STRING_id == x],
                                        USE.NAMES = FALSE)
  
resp_string2$RESP <- resp_string2$Gene %in% resp$Gene
resp_string2 <- left_join(resp_string2, resp[,c("id", "Gene", "shift")] %>%
                           dplyr::rename(shift_RESP = shift),
                         by = c("id", "Gene"))


## Plotting full graph if necessary
colnames(resp_string2)[1] <- "name"
string_graph <- igraph::graph_from_data_frame(resp_string_interactions2, 
                                              vertices = resp_string2[which(!is.na(match(resp_string2$name,
                                                                                        unique(c(resp_string_interactions2$from,
                                                                                                 resp_string_interactions2$to))
                                                                                        ))),],
                                              directed = FALSE
                                              ) 
set.seed(123)
p <- ggraph(string_graph, layout = 'fr', weights = combined_score,
            niter = 10000) + 
  geom_edge_link(linewidth = 0.8, aes(width = combined_score, alpha = combined_score), 
                 show.legend = FALSE) +
  geom_node_point(aes(shape = RESP, color = QP_Venetoclax12h, fill = IS_Ventoclax12h), 
                  size = 7, stroke = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        name = "QP 12h", na.value = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       name = "I-score 12h", na.value = "black") + 
  new_scale_fill() + new_scale_color() +
  geom_node_point(aes(shape = RESP, color = shift_RESP, fill = shift_RESP, filter = !is.na(shift_RESP)), 
                 size = 7, stroke = 2) +
  scale_color_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c"),
                     guide = "none", na.value = "black") +
  scale_fill_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c"),
                    name = "RESP shift", na.value = "black",
                    guide = guide_legend(override.aes = list(shape = 23, stroke = 0, size = 7))) +
  geom_node_label(aes(x = x, y = y, label = name), 
                  max.time = 2, max.iter = 20000,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  #angle = 0, hjust = 0.5, nudge_y = 0.1, 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 10) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), name = "", 
                     labels = c("interactors", "RESP"),
                     guide = guide_legend(override.aes = list(size = 5))) +
  labs(title = "STRING network from RESP proteins of the Astral MS Venetoclax dataset") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 20),
        plot.subtitle = element_text(size = 18))

#ggsave("Figures/STRING_network_RESP_AstralVenetoclax - 12h.png", 
#       plot = p,  dpi = 300, device = "png", width = 26, height = 16)


### graph clustering
resp_string_graph_wght <- resp_string_interactions2$combined_score
resp_string_clusters <- igraph::cluster_walktrap(string_graph, steps = 5,
                                                 weights = resp_string_graph_wght)

resp_string_clusters <- sapply(1:max(resp_string_clusters$membership), 
                               function(com) resp_string_clusters$names[which(resp_string_clusters$membership == com)])
length(resp_string_clusters) # 15 clusters
sum(sapply(resp_string_clusters, length) > 2 & 
      sapply(resp_string_clusters, function(z) any(z %in% resp$Gene))) 
# 9 clusters will be kept

## Using GSAI to concisely name clusters (https://idekerlab.ucsd.edu/gsai/)
## --> using GPT4o and their prompt; making the same query to GPT, get back results
## --> GSAI from R 
GSAI_naming <- function(proteins){
  gsai_api <- "sk-hnkJlW4jWpHLmymDYeoJbJPVS77QuogiUVOMhBKcWjT3BlbkFJyHt-Jd4mzrz192vMrlNXJ5rQ-QHhkWqTiYs8VjpdsA"
  gsai_prompt <- 
    "Write a critical analysis of the biological processes performed by this system of interacting proteins.
 Base your analysis on prior knowledge available in your training data. After completing your analysis, propose a brief and detailed name for the most prominent biological process performed by the system.
 After completing your analysis, please also assign a confidence level to the process name you selected. This confidence level should follow the name in parentheses and be one of the following: Low, Medium, or High. A low confidence level indicates the least confidence, while a high confidence level reflects the most confidence. This confidence level helps gauge how accurately the chosen name represents the functions and activities within the system of interacting proteins. When determining your confidence level, consider the proportion of genes in the protein system that participate in the identified biological process. For instance, if you select “Ribosome biogenesis” as the process name but only a few genes in the system contribute to this process, the confidence level should be low compared to a scenario where a majority of the genes are involved in “Ribosome biogenesis.”
 Put your chosen name at the top of the analysis as 'Name: <name>'. 
 Put your confidence in a new line after the name as ‘LLM self-assessed confidence: <confidence level>’
 Be concise: Avoid unnecessary words.
 Be factual: Do not editorialize.
 Be specific: Avoid overly general statements such as ‘the proteins are involved in various cellular processes’.
 Group proteins: If group of  proteins has similar functions then discuss their interplay, synergistic, or antagonistic effects, and functional integration within the system, instead of listing facts about individual proteins
 Avoid generic process names: Choose detailed and specific names, avoiding generic names like ‘Cellular Signaling and Regulation’.
 Contextual relevance: Ensure the analysis is relevant to the specific biological context of the system of interacting proteins if provided.

 If you cannot identify a prominent biological process for the proteins in the system, I want you to communicate this in you analysis and name the process: 'System of unrelated proteins'. Provide a confidence of ‘None’ for a 'System of unrelated proteins'.
    
 To help you in your work, I am providing an example system of interacting proteins and the corresponding example analysis output.

 The example system of interacting proteins is:
 PDX1, SLC2A2, NKX6-1, GLP1, GCG.

  The example analysis output is:
  Name: Pancreatic development and glucose homeostasis 
  LLM self-assessed confidence: High
  
  1. PDX1 is a homeodomain transcription factor involved in the specification of the early pancreatic epithelium and 
  its subsequent differentiation. It activates the transcription of several genes including insulin, somatostatin, glucokinase 
  and glucose transporter type 2. It is essential for maintenance of the normal hormone-producing phenotype in the 
  pancreatic beta-cell. In pancreatic acinar cells, forms a complex with PBX1b and MEIS2b and mediates the activation 
  of the ELA1 enhancer.
  
  2. NKX6-1 is also a transcription factor involved in the development of pancreatic beta-cells during the secondary transition. 
  Together with NKX2-2 and IRX3, controls the generation of motor neurons in the neural tube and belongs to the neural progenitor 
  factors induced by Sonic Hedgehog (SHH) signals.
  
  3.GCG and GLP1, respectively glucagon and glucagon-like peptide 1, are involved in glucose metabolism and homeostasis. 
  GCG raises blood glucose levels by promoting gluconeogenesis and is the counter regulatory hormone of Insulin. 
  GLP1 is a potent stimulator of Glucose-Induced Insulin Secretion (GSIS). Plays roles in gastric motility and 
  suppresses blood glucagon levels. Promotes growth of the intestinal epithelium and pancreatic islet mass both by islet 
  neogenesis and islet cell proliferation.
  
  4. SLC2A2, also known as GLUT2, is a facilitative hexose transporter. In hepatocytes, it mediates bi-directional 
  transport of glucose accross the plasma membranes, while in the pancreatic beta-cell, it is the main transporter responsible 
  for glucose uptake and part of the cell's glucose-sensing mechanism. It is involved in glucose transport in the small intestine 
  and kidney too.
  
  To summarize, the genes in this set are involved in the specification, differentiation, growth and functionality of the pancreas, 
  with a particular emphasis on the pancreatic beta-cell. Particularly, the architecture of the pancreatic islet ensures proper 
  glucose sensing and homeostasis via a number of different hormones and receptors that can elicit both synergistic and antagonistic 
  effects in the pancreas itself and other peripheral tissues.
      
  
  Here are the interacting proteins:

  Proteins:
  "
  
  protein_list <- paste(proteins, collapse = ",")
  
  gsai_response <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    httr::add_headers(
      Authorization = paste("Bearer", gsai_api)
    ),
    httr::content_type_json(),
    encode = "json",
    body = list(
      model = "gpt-4o-2024-08-06",
      messages = list(
        list(role = "system", content = "Process Naming and Analysis"),
        list(role = "user", content = paste(gsai_prompt, protein_list))
      )
    )
  )
  
  gsai_response <- httr::content(gsai_response, "text", encoding = "UTF-8")
  gsai_response <- sub('.*content\\": \\"Name: ', "", gsai_response)
  
  gsai_name <- gsub(" {1,}$", "", sub("\\\\.*", "", gsai_response))
  gsai_confidence <- sub("\\\\.*", "", sub(".*LLM self-assessed confidence: ", "", gsai_response))
  
  return(data.frame(GSAI_name = gsai_name, GSAI_confidence = gsai_confidence))
}
# apply the function to each cluster to name them
resp_string_clusters_GSAIwconf <- lapply(resp_string_clusters, GSAI_naming)
resp_string_clusters_GSAI <- lapply(resp_string_clusters_GSAIwconf, "[[", 1)
resp_string_clusters_GSAI <- lapply(resp_string_clusters_GSAI, gsub, pattern = ", |, and ",  replacement = ",\n")
resp_string_clusters_GSAI <- lapply(resp_string_clusters_GSAI, gsub, pattern = " and ", replacement = "\nand ")

### !! Name might slightly differ between run !! (more consistent with bigger clusters) ###

resp_string_all <- list()
resp_string_interactions_all <- list()
for(string_cluster in 1:length(resp_string_clusters)){
  resp_string2 <- resp_string[which(!is.na(match(resp_string$Gene, resp_string_clusters[[string_cluster]]))),]
  resp_string_interactions2 <- resp_string_interactions[which(!is.na(match(resp_string_interactions$from, resp_string2$STRING_id)) &
                                                                !is.na(match(resp_string_interactions$to, resp_string2$STRING_id))
                                                              ),]
  
  resp_string_interactions2$from <- sapply(resp_string_interactions2$from,
                                           function(x) resp_string2$Gene[resp_string2$STRING_id == x],
                                           USE.NAMES = FALSE)
  resp_string_interactions2$to <- sapply(resp_string_interactions2$to,
                                         function(x) resp_string2$Gene[resp_string2$STRING_id == x],
                                         USE.NAMES = FALSE)
  
  resp_string2$RESP <- resp_string2$Gene %in% resp$Gene
  
  resp_string2 <- left_join(resp_string2, resp[,c("id", "Gene", "shift")] %>%
                              dplyr::rename(shift_RESP = shift),
                            by = c("id", "Gene"))
  colnames(resp_string2)[1] <- "name"
  
  
  resp_string2$cluster <- paste0("Cluster ", string_cluster, "\n",
                                 resp_string_clusters_GSAI[[string_cluster]])
  resp_string2$cluster_confidence_naming <- resp_string_clusters_GSAIwconf[[string_cluster]]$GSAI_confidence
  
  if(nrow(resp_string_interactions2)){
    resp_string_interactions2$cluster <-  paste0("Cluster ", string_cluster, "\n",
                                                 resp_string_clusters_GSAI[[string_cluster]])
    resp_string_interactions2$cluster_confidence_naming <- resp_string_clusters_GSAIwconf[[string_cluster]]$GSAI_confidence
  }
  
  resp_string_all[[string_cluster]] <- resp_string2
  resp_string_interactions_all[[string_cluster]] <- resp_string_interactions2
}


### Plotting graphs
vertices <- as.data.frame(Reduce(rbind, resp_string_all))
relations <- as.clusterrelations <- as.data.frame(Reduce(rbind, resp_string_interactions_all))

## saving in excel file
openxlsx::write.xlsx(list("RESP_interactors_GSAI" = vertices,
                          "RESP_interactors_relations" = relations),
                     "Results/VenetoclaxAstral_RESP_STRINGdbInteractorsGSAI-Fig3.xlsx")

## filter clusters
vertices <- vertices[which(!is.na(match(vertices$cluster,
                                        names(which(table(vertices$cluster) > 2))
                                        ))),]
# only keeping clusters with RESP proteins:
vertices %>% group_by(cluster) %>% summarise(anyRESP = any(RESP)) %>% dplyr::select(anyRESP) %>% all()
# TRUE

relations <- relations[which(!is.na(match(relations$cluster,
                                          unique(vertices$cluster)
                                          ))),]

string_graph2 <- igraph::graph_from_data_frame(relations, 
                                               vertices = vertices,
                                               directed = FALSE
                                               ) %>%
  igraph::set_vertex_attr("cluster", value = factor(vertices$cluster, c(unique(vertices$cluster)[c(1,3,6)],
                                                                        unique(vertices$cluster)[c(4,5,2)],
                                                                        unique(vertices$cluster)[c(7,9,8)]
                                                                        )
                                                    )
                          ) %>%
  igraph::set_edge_attr("cluster", value = factor(relations$cluster, c(unique(vertices$cluster)[c(1,3,6)],
                                                                       unique(vertices$cluster)[c(4,5,2)],
                                                                       unique(vertices$cluster)[c(7,9,8)]
                                                                       )
                                                  )
                        )

set.seed(123)
p <- ggraph(string_graph2, layout = 'fr', weights = combined_score,
            niter = 10000) + 
  geom_edge_link(linewidth = 0.8, aes(width = combined_score, alpha = combined_score), 
                 show.legend = FALSE) +
  geom_node_point(aes(shape = RESP, color = QP_Venetoclax12h, fill = IS_Ventoclax12h), 
                  size = 8, stroke = 2.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        name = "QP 12h (border)", na.value = "black",
                        guide = guide_colorbar(order = 4, barheight = unit(3, "inches"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       name = "I-score 12h (inside)", na.value = "black",
                       guide = guide_colorbar(order = 3, barheight = unit(3, "inches"))) + 
  new_scale_fill() + new_scale_color() +
  geom_node_point(aes(shape = RESP, color = shift_RESP, fill = shift_RESP, filter = !is.na(shift_RESP)), 
                  size = 8, stroke = 2.5) +
  scale_color_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c"),
                     guide = "none", na.value = "black") +
  scale_fill_manual(values = c("3h" = "#ff8080", "12h" = "#d090f8", "3h & 12h"  = "#e7487c"),
                    name = "RESP shift", na.value = "black",
                    guide = guide_legend(override.aes = list(shape = 23, stroke = 0, size = 9), order = 2)) +
  geom_node_label(aes(x = x, y = y, label = name, filter = !RESP), 
                  max.time = 2, max.iter = 20000, size = 5,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 3) +
  geom_node_label(aes(x = x, y = y, label = name, filter = RESP), 
                  max.time = 2, max.iter = 20000, size = 7,
                  color = "black", show.legend = FALSE, fontface = "bold", 
                  repel = TRUE, min.segment.length = 0, max.overlaps = 20) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), name = "", 
                     labels = c("CETSA/QP", "RESP"), 
                     guide = guide_legend(override.aes = list(size = 9), order = 1)) +
  facet_wrap(~cluster, scales = "free") +
  labs(title = "Clustered STRING networks of the RESP proteins\nand their top 50 interactors shifting in CETSA/QP\nof the Astral MS Venetoclax dataset") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 30),
        panel.border = element_rect(color = "black"),
        strip.background = element_rect(fill = "white", color = "black"), 
        strip.text = element_text(face = "bold", color = "black", size = 26))

ggsave("Figures/Fig 3 STRING_network_RESP_AstralVenetoclax - all clusters - 12h.png",
       plot = p,  dpi = 300, device = "png", width = 32, height = 28)

## data cleaning
rm(resp, prot, hits, string_db, string_db_proteins, string_idbase, nmin_interactors, string_id, 
   resp_string, resp_string_interactions, resp_string2, resp_string_interactions2, string_graph, 
   resp_string_graph_wght, resp_string_clusters, resp_string_clusters_GSAIwconf, resp_string_clusters_GSAI, 
   resp_string_all, resp_string_interactions_all, string_cluster, vertices, relations, string_graph2, p)
