library(IMPRINTS.CETSA.app)


# review time categorization of each RESP
filteredresp_timelist <- list("RESP_3h" = unique(vene_astral_filteredresp$id[vene_astral_filteredresp$treatment == "Venetoclax3h"]),
                              "RESP_12h" = unique(vene_astral_filteredresp$id[vene_astral_filteredresp$treatment == "Venetoclax12h"])
                              )

filteredresp_timelist <- com_protein_loop(filteredresp_timelist)
sapply(filteredresp_timelist, length)
# RESP_3h  RESP_12h RESP_3h & RESP_12h 
#      49       27                 19 

# several reason can cause a non overlap:
# data is slightly more noisy in another time point and categorization didn't give the same result
# p-value and RESP-score distribution being different, some proteins can fall just below the cutoff but can still show a RESP effect
# cleavage indeed occur later or earlier or just too noisy to conclude anything

### let's check the 49 proteins with a RESP effect at 3h at 12h
### and th 27 at 12h at 3h

## Different categorization
## if categorized as other than Unconfident, SP or MP can occur to degradation so they will be counted as RESP
## if Unconfident need to be checked manually
# RESP 3h different category at 12h:
resp3h_catdiff12h <- vene_astral_resp12h[which(!is.na(match(vene_astral_resp12h$id,
                                       filteredresp_timelist$RESP_3h[filteredresp_timelist$RESP_3h %in% 
                                                                       vene_astral_resp12h$id])
                                 )),]
resp3h_catdiff12h %>%
  group_by(category) %>%
  summarise(Gene = paste(Gene, collapse = " "))
# MP --> SUGP1      
# SP --> DKC1 DDX21 PLEC LARP1 
# Unconfident --> SF3B1 NOLC1 PRPF40A RBBP6 STRN HSPA4L GAK
# SF3B1 | most likely start to be degraded (several C-term side peptides with very low FC) --> RESP
# NOLC1 | most likely due to PTM
# PRPF40A | C-term side more noisy at 12h --> RESP
# RBBP6 | caused by 3 N-term peptides with slightly greater FC but still clear RESP --> RESP
# STRN | most likely start to be degraded (several N-term side peptides with very low FC) --> RESP
# HSPA4L | should be RESP, unconfident most likely due to 622-632 --> RESP
# GAK | very similar profile at 3h and 12h, N-term might start to be degraded --> RESP

# Only NOLC1 stays only 3h
# SUGP1, DKC1, DDX21, PLEC, LARP1, SF3B1, PRPF40A, RBBP6, STRN, HSPA4L, GAK pass to 3h and 12h
resp3h_recat_to3hand12h <- resp3h_catdiff12h$id[which(!is.na(match(resp3h_catdiff12h$Gene,
                                                              c("SUGP1", "DKC1", "DDX21", "PLEC", "LARP1", 
                                                                "SF3B1", "PRPF40A", "RBBP6", "STRN", 
                                                                "HSPA4L", "GAK")
                                                              ))
                                                 )]

## RESP 12h different category at 3h:
resp12h_catdiff3h <- vene_astral_resp3h[which(!is.na(match(vene_astral_resp3h$id,
                                       filteredresp_timelist$RESP_12h[filteredresp_timelist$RESP_12h %in% 
                                                                       vene_astral_resp3h$id])
                                 )),]
resp12h_catdiff3h %>%
  group_by(category) %>%
  summarise(Gene = paste(Gene, collapse = " "))
# Unconfident --> HNRNPA3 and CCNT1
# HNRNPA3 | cannot be considered a RESP at 3h, 21-29 shift at 3h and doesn't at 12h
# CCNT1 | less prominent at 3h than at 12h but would consider it also a RESP at 3h

# HNRNPA3 stays only at 12h 
# CCNT1 pass to 3h and 12h

resp12h_recat_to3hand12h <- resp12h_catdiff3h$id[which(!is.na(match(resp12h_catdiff3h$Gene,
                                                              c("CCNT1")
                                                              ))
                                                 )]

## Non hits
## RESP 3h that didn't pass cutoffs at 12h
resp3h_non12h <- openxlsx::read.xlsx(list.files(path = list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_venetoclax_astral_12h"),
                                                pattern = "full", full.names = TRUE)
                                     ) %>%
  dplyr::filter(id %in% filteredresp_timelist$RESP_3h & 
                  !(id %in% vene_astral_resp12h$id))
# recategorize
resp3h_non12h <- imprints_categorize_peptides(vene_astral_pepnorm12h, resp3h_non12h, 
                                              "Vehicle12h", FALSE)
# plot
imprints_barplotting_categorize_peptides(vene_astral_pepnorm12h, resp3h_non12h,
                                         "Venetoclax12h", "Vehicle12h", 
                                         "peptide_one", "purple",
                                         "VenetoclaxAstral_RESP3h_non12h_barplotting")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP3h_non12h_barplotting.pdf"),
            "Results/VenetoclaxAstral_RESP3h_non12h_barplotting.pdf")

# if same cleavage site and p-value below 0.05, will consider also a RESP at 12h
resp3h_recat_to3hand12h <- c(resp3h_recat_to3hand12h,
                             resp3h_non12h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                               left_join(vene_astral_resp3h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                         by = c("id", "Gene"), suffix = c(".12h", ".3h")) %>%
                               dplyr::filter(cleaved_site.12h == cleaved_site.3h &
                                               combined_pvalue <= 0.05) %>%
                               dplyr::select(id) %>% c() %>% unlist() %>% unname()
                             )

# if different cleavage site and p-value below 0.05
# will consider also a RESP at 12h if categorized as RESP with RESP score > 0.2 
# or if not categorized as RESP, required RESP_score > 0.8
resp3h_recat_to3hand12h <- c(resp3h_recat_to3hand12h,
                             resp3h_non12h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                               left_join(vene_astral_resp3h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                         by = c("id", "Gene"), suffix = c(".12h", ".3h")) %>%
                               dplyr::filter(cleaved_site.12h != cleaved_site.3h & combined_pvalue <= 0.05 & 
                                               (category == "RESP" & abs(RESP_score.12h) >= 0.2 | 
                                                  abs(RESP_score.12h) >= 0.8)
                                             )  %>%
                               dplyr::select(id) %>% c() %>% unlist() %>% unname()
                             )

# if different cleavage site and p-value greater 0.05
# required RESP_score > 1
resp3h_recat_to3hand12h <- c(resp3h_recat_to3hand12h,
                             resp3h_non12h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                               left_join(vene_astral_resp3h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                         by = c("id", "Gene"), suffix = c(".12h", ".3h")) %>%
                               dplyr::filter(cleaved_site.12h != cleaved_site.3h & combined_pvalue > 0.05 &
                                               abs(RESP_score.12h) >= 1)  %>%
                               dplyr::select(id) %>% c() %>% unlist() %>% unname()
                             )

## RESP 12h that didn't pass cutoffs at 3h
resp12h_non3h <- openxlsx::read.xlsx(list.files(path = list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_venetoclax_astral_3h"),
                                                pattern = "full", full.names = TRUE)) %>%
  dplyr::filter(id %in% filteredresp_timelist$RESP_12h & 
                  !(id %in% vene_astral_resp3h$id))

# recategorize
resp12h_non3h <- imprints_categorize_peptides(vene_astral_pepnorm3h, resp12h_non3h, 
                                              "Vehicle3h", FALSE)
# plot
imprints_barplotting_categorize_peptides(vene_astral_pepnorm3h, resp12h_non3h,
                                         "Venetoclax3h", "Vehicle3h", 
                                         "peptide_one", "red",
                                         "VenetoclaxAstral_RESP12h_non3h_barplotting")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP12h_non3h_barplotting.pdf"),
            "Results/VenetoclaxAstral_RESP12h_non3h_barplotting.pdf")

# if same cleavage site and p-value below 0.05, will consider also a RESP at 3h
resp12h_recat_to3hand12h <- c(resp12h_recat_to3hand12h, 
                              resp12h_non3h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                                left_join(vene_astral_resp12h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                          by = c("id", "Gene"), suffix = c(".3h", ".12h")) %>%
                                dplyr::filter(cleaved_site.3h == cleaved_site.12h &
                                                combined_pvalue <= 0.05)  %>%
                                dplyr::select(id) %>% c() %>% unlist() %>% unname()
                              )

# if different cleavage site and p-value below 0.05
# will consider also a RESP at 3h if categorized as RESP with RESP score > 0.2 
# or if not categorized as RESP, required RESP_score > 0.8
resp12h_recat_to3hand12h <- c(resp12h_recat_to3hand12h,
                              resp12h_non3h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                               left_join(vene_astral_resp12h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                         by = c("id", "Gene"), suffix = c(".3h", ".12h")) %>%
                               dplyr::filter(cleaved_site.12h != cleaved_site.3h & combined_pvalue <= 0.05 & 
                                               (category == "RESP" & abs(RESP_score.3h) >= 0.2 | 
                                                  abs(RESP_score.3h) >= 0.8)
                                             )  %>%
                               dplyr::select(id) %>% c() %>% unlist() %>% unname()
                              )

# if different cleavage site and p-value greater 0.05
# required RESP_score > 1
resp12h_recat_to3hand12h <- c(resp12h_recat_to3hand12h,
                              resp12h_non3h[,c("id", "Gene", "cleaved_site", "combined_pvalue", "RESP_score", "category")] %>%
                               left_join(vene_astral_resp12h[,c("id", "Gene", "cleaved_site", "RESP_score")],
                                         by = c("id", "Gene"), suffix = c(".3h", ".12h")) %>%
                               dplyr::filter(cleaved_site.12h != cleaved_site.3h & combined_pvalue > 0.05 &
                                               abs(RESP_score.3h) >= 1)  %>%
                               dplyr::select(id) %>% c() %>% unlist() %>% unname()
                              )

### reupdating list
filteredresp_timelist <- list("RESP_3h" = c(unique(vene_astral_filteredresp$id[vene_astral_filteredresp$treatment == "Venetoclax3h"]),
                                            resp12h_recat_to3hand12h),
                              "RESP_12h" = c(unique(vene_astral_filteredresp$id[vene_astral_filteredresp$treatment == "Venetoclax12h"]),
                                             resp3h_recat_to3hand12h)
                              )
sapply(com_protein_loop(filteredresp_timelist), length)
# RESP_3h   RESP_12h RESP_3h & RESP_12h 
#      14        11                 70 

## Final Venn diagram         
names(filteredresp_timelist) <- c("3h", "12h")
VennDiagram::venn.diagram(filteredresp_timelist, 
                          "Figures/VennDiagram_VenetoclaxAstral_catasRESP_3h12h_recategorized.png",
                          output = TRUE, disable.logging = TRUE, imagetype = "png",
                          fill = c("red", "purple"),
                          fontface = "bold", margin = 0.075,
                          fontfamiliy = "sans",
                          cat.cex = 1.8, cex = 1.7, lwd = 2,
                          cat.fontface = "bold",
                          cat.fontfamily = "sans")

### adding all 95 proteins in one table
vene_astral_resp_final <- Reduce(rbind,
       list(vene_astral_filteredresp %>%
              mutate(significant_shift = TRUE),
            resp3h_catdiff12h[,colnames(vene_astral_filteredresp)] %>%
              mutate(significant_shift = TRUE),
            resp12h_catdiff3h[,colnames(vene_astral_filteredresp)] %>%
              mutate(significant_shift = TRUE),
            resp3h_non12h[,colnames(vene_astral_filteredresp)] %>%
              mutate(significant_shift = FALSE),
            resp12h_non3h[,colnames(vene_astral_filteredresp)] %>%
              mutate(significant_shift = FALSE)
            )
       )
vene_astral_resp_final$time_category <- sapply(vene_astral_resp_final$id, 
                                               function(p)
                                                 paste(names(which(sapply(filteredresp_timelist, 
                                                                          function(x) p %in% x))), 
                                                       collapse = " & "),
                                               USE.NAMES = FALSE
                                               )
vene_astral_resp_final <- vene_astral_resp_final[order(vene_astral_resp_final$Gene),]


