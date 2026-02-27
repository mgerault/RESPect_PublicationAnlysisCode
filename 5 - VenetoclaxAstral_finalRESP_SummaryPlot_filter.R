library(IMPRINTS.CETSA.app)

## modifying plotting function
barplotting_peptides <- function(x, color){
  require(ggtext)
  
  # getting cleavage info
  x_cleavageinfo <- vene_astral_resp_final[which(vene_astral_resp_final$id == 
                                                   x$Master.Protein.Accessions[1]),] 
  
  x_cleavageinfo <- paste0("<br>**At 3h**: RESP-score=", 
         round(unique(x_cleavageinfo$RESP_score[which(x_cleavageinfo$treatment == "Venetoclax3h")]), 3),
         " p-value=",
         signif(unique(x_cleavageinfo$combined_pvalue[which(x_cleavageinfo$treatment == "Venetoclax3h")]), 3),
         ifelse(unique(x_cleavageinfo$significant_shift[which(x_cleavageinfo$treatment == "Venetoclax3h")]),
                " with", " without"),
         " a significant shift and categorized as ",
         unique(x_cleavageinfo$category[which(x_cleavageinfo$treatment == "Venetoclax3h")]),
         
         "<br>**At 12h**: RESP-score=", 
         round(unique(x_cleavageinfo$RESP_score[which(x_cleavageinfo$treatment == "Venetoclax12h")]), 3),
         " p-value=",
         signif(unique(x_cleavageinfo$combined_pvalue[which(x_cleavageinfo$treatment == "Venetoclax12h")]), 3),
         ifelse(unique(x_cleavageinfo$significant_shift[which(x_cleavageinfo$treatment == "Venetoclax12h")]),
                " with", " without"),
         " a significant shift and categorized as ",
         unique(x_cleavageinfo$category[which(x_cleavageinfo$treatment == "Venetoclax12h")]),
         
         "<br><br>**", unique(x_cleavageinfo$Gene), "** was categorized as having a RESP at ", unique(x_cleavageinfo$time_category),
         ifelse(length(unique(x_cleavageinfo$cleaved_site)) == 1,
                paste0(" with the same computed cleavage site, **", unique(x_cleavageinfo$cleaved_site),
                       "**, at 3 and 12 hours resulting in the following sequence (+/- 10AA):<br>",
                       gsub("D", "<b style='color:red'>D</b>", unique(x_cleavageinfo$cleaved_sequence))
                ),
                paste0(" with the computed cleavage sites **", 
                       unique(x_cleavageinfo$cleaved_site[which(x_cleavageinfo$treatment == "Venetoclax3h")]),
                       "** at 3h and **",
                       unique(x_cleavageinfo$cleaved_site[which(x_cleavageinfo$treatment == "Venetoclax12h")]),
                       "** at 12h  resulting in the following sequences (+/- 10AA):<br>",
                       gsub("D", "<b style='color:red'>D</b>",
                            unique(x_cleavageinfo$cleaved_sequence[which(x_cleavageinfo$treatment == "Venetoclax3h")])
                       ),
                       " & ",
                       gsub("D", "<b style='color:red'>D</b>",
                            unique(x_cleavageinfo$cleaved_sequence[which(x_cleavageinfo$treatment == "Venetoclax12h")])
                       )
                )
         ),
         "<br><br>**", unique(x_cleavageinfo$Gene), "** was",
         ifelse(any(x_cleavageinfo$in_merCASBA),
                paste0(" reported to be cleaved by merCASBA at **",
                       paste(unique(paste0(sub(".{3}", "", 
                                               gsub("[a-z]", "", 
                                                    x_cleavageinfo$cleaved_sequence.merCASBA)), 
                                           x_cleavageinfo$cleaved_site.merCASBA)), 
                             collapse = " & "),
                       "**, reference(s): ", unique(x_cleavageinfo$merCASBA_reference),
                       ifelse(any(!is.na(x_cleavageinfo$merCASBA_PMID)),
                              paste0(" (PMID: ",
                                     unique(na.omit(x_cleavageinfo$merCASBA_PMID)), ")"),
                              ""),
                       "<br>The cleavage site does ",
                       ifelse(any(x_cleavageinfo$cleaved_site_in_merCASBA, na.rm = TRUE),
                              paste0("correspond to the one computed at ",
                                     paste(gsub("Venetoclax", "",
                                                unique(x_cleavageinfo$treatment[which(x_cleavageinfo$cleaved_site_in_merCASBA)])
                                     ),
                                     collapse = " & ")
                              ),
                              "not correspond to the one computed at 3h & 12h"
                       )
                ),
                " not reported to be cleaved by merCASBA"),
         "<br>"
         )
  x$countNum <- NULL
  
  # preparing data - mean and se
  x <- x %>%
    tidyr::gather("key", "value", -Master.Protein.Accessions, -description,
                  -Positions.in.Master.Proteins,
                  -Annotated.Sequence, -Modifications) %>%
    tidyr::separate(key, into = c("temperature", "rep", "treatment"), sep = "_") %>%
    group_by(Master.Protein.Accessions, description,
             Positions.in.Master.Proteins,
             Annotated.Sequence, Modifications,
             treatment, temperature) %>%
    summarise(se = sd(value, na.rm = TRUE)/sqrt(length(na.omit(value))),
              value = mean(value, na.rm = TRUE))
  
  # preparing title
  x$Gene <- sapply(x$description, IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE)
  x$description <- sapply(x$description, IMPRINTS.CETSA.app:::getProteinName, USE.NAMES = FALSE)
  x$temperature <- factor(x$temperature)
  x$Positions.in.Master.Proteins <- sub(".* ", "", x$Positions.in.Master.Proteins)
  x$Positions.in.Master.Proteins <- gsub("\\[|\\]", "", x$Positions.in.Master.Proteins)
  
  # formating modification
  x$Modifications <- unlist(lapply(strsplit(x$Modifications, "\\];"),
                                   function(y){
                                     y <- y[-grep("TMT", y)]
                                     if(length(y)){
                                       y <- sub("\\d{1}x", "", y)
                                       y <- sub("^ ", "", y)
                                       y <- gsub("\\[|\\]", "", y)
                                       y <- paste(y, collapse = "\n")
                                     }
                                     else{
                                       y <- ""
                                     };
                                     y
                                   })
  )
  
  # formating labels/title
  x$Positions.in.Master.Proteins <- paste0(x$Positions.in.Master.Proteins, "\n",
                                           x$Modifications)
  x$Positions.in.Master.Proteins <- factor(x$Positions.in.Master.Proteins,
                                           unique(x$Positions.in.Master.Proteins)[order(sapply(strsplit(gsub("\\[|\\]|\n.*", "",
                                                                                                             unique(x$Positions.in.Master.Proteins)
                                           ),
                                           "-"),
                                           function(y) sum(as.numeric(y)))
                                           )])
  
  # handling if QP in temperature
  x$temperature <- factor(x$temperature)
  x$QP <- FALSE
  lvl_temperature <- levels(x$temperature)
  if("36C" %in% x$temperature){
    x$QP[which(x$temperature == "36C")] <- TRUE
    x$temperature <- factor(sub("36C", "QP", x$temperature),
                            levels = c("QP", lvl_temperature[-grep("36C", lvl_temperature)]))
  }
  
  x$treatment <- factor(x$treatment, levels = c("Venetoclax3h", "Venetoclax12h"))
  x$condition <- paste0(x$temperature, "_", x$treatment)
  x$condition <- factor(x$condition, 
                        levels = c(paste0(levels(x$temperature), "_", levels(x$treatment)[1]),
                                   paste0(levels(x$temperature), "_", levels(x$treatment)[2])))
  
  # plot
  g <- ggplot(x, aes(condition, value, fill = treatment)) +
    geom_bar(stat = "identity", aes(color = QP)) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.1) +
    facet_wrap(~Positions.in.Master.Proteins, strip.position = "bottom",
               nrow = 1) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
    scale_x_discrete(labels = gsub("_.{1,}", "", levels(x$condition))) + 
    guides(color = "none") +
    labs(y = "fold-change (log2)", x = "peptides",
         title = paste(x$Master.Protein.Accessions[1], "\n",
                       x$description[1], "\n",
                       x$Gene[1]),
         subtitle = x_cleavageinfo) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1.75), face = "bold"),
          plot.subtitle = element_markdown(size = rel(1.5)),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(face = "bold", size = rel(1.4)),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.9)),
          axis.text.y = element_text(size = rel(1.8)),
          legend.text = element_text(face = "bold", size = rel(1.4)),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold", angle = 90, size = rel(1.4)))
  
  return(g)
}

# join datasets to plot both time point
vene_astral_resp_final_plt <- full_join(vene_astral_pepdiff3h[which(!is.na(match(vene_astral_pepdiff3h$Master.Protein.Accessions,
                                                                                 unique(vene_astral_resp_final$id)
                                                                                 )
                                                                           )),
                                                              -grep("_Vehicle3h$|countNum", colnames(vene_astral_pepdiff3h))],
                                        vene_astral_pepdiff12h[which(!is.na(match(vene_astral_pepdiff12h$Master.Protein.Accessions,
                                                                                  unique(vene_astral_resp_final$id)
                                                                                  )
                                                                            )),
                                                               -grep("_Vehicle12h$", colnames(vene_astral_pepdiff12h))],
                                        by = colnames(vene_astral_pepdiff3h)[1:5])

# apply function to all protein of interest i.e. the filtered RESP
l <- plyr::dlply(vene_astral_resp_final_plt, 
                 plyr::.(Master.Protein.Accessions), 
                 .fun = barplotting_peptides, 
                 color = c("red", "purple"))

# order plot accoring to Genes
l <- l[order(sapply(names(l), 
                    function(x) unique(vene_astral_resp_final$Gene[which(vene_astral_resp_final$id == x)]),
                    USE.NAMES = FALSE))]


# save in a pdf
class(l) <- c("arrangelist", "ggplot", class(l))
ggpubr::ggexport(filename = "Results/VenetoclaxAstral_RESPfinal_SummaryPlot.pdf",
                 plotlist = l, height = 12, width = 42)

rm(l, vene_astral_resp_final_plt)


### Removing potential false positive after manual inspection
vene_astral_resp_potentialFP <- unique(vene_astral_resp_final$id[which(!is.na(match(vene_astral_resp_final$Gene,
                                                                                    c("KIF21A", "NOP14", "USP36")
                                                                                    )
                                                                              ))])
# careful with: PRDX4, BAG6, RBM27, NAT10 and MTA1 as only two peptides driving the change
# --> Often too small regions for giving folded domain, I.e no hydrophobic core, I.e no CETSA shift.
# but only PRDX4 and BAG6 have a region of only ~30/40 AA but BAG6 and PRDX4 do seem to have an alpha-helix at these regions
# will keep them all
filteredresp_timelist_refilter <- list("3h" = setdiff(unique(vene_astral_resp_final$id[grep("3h", vene_astral_resp_final$time_category)]),
                                                      vene_astral_resp_potentialFP),
                                       "12h" = setdiff(unique(vene_astral_resp_final$id[grep("12h", vene_astral_resp_final$time_category)]),
                                                       vene_astral_resp_potentialFP)
                                       )
sapply(com_protein_loop(filteredresp_timelist_refilter), length)
# RESP_3h   RESP_12h RESP_3h & RESP_12h 
#      14        9                 69 
## --> 92 RESP

## Final Venn diagram         
VennDiagram::venn.diagram(filteredresp_timelist_refilter, 
                          "Figure/VennDiagram_VenetoclaxAstral_catasRESP_3h12h_recategorized_refiltered.png",
                          output = TRUE, disable.logging = TRUE, imagetype = "png",
                          fill = c("red", "purple"),
                          fontface = "bold", margin = 0.075,
                          fontfamiliy = "sans",
                          cat.cex = 1.8, cex = 1.7, lwd = 2,
                          cat.fontface = "bold",
                          cat.fontfamily = "sans")

length(unique(vene_astral_resp_final$id)) - length(vene_astral_resp_potentialFP) # 92 unique proteins
vene_astral_resp_final %>% 
  dplyr::filter(!(id %in% vene_astral_resp_potentialFP)) %>%
  group_by(id, Gene) %>%
  summarise(in_merCASBA = any(in_merCASBA),
            site_in_merCASBA = any(cleaved_site_in_merCASBA)) %>%
  ungroup() %>%
  summarise(in_merCASBA = sum(in_merCASBA), # 36 / 92
            site_in_merCASBA = sum(site_in_merCASBA) # 9 / 36
  )

vene_astral_resp_final$kept_for_further_analysis <- !(vene_astral_resp_final$id %in% vene_astral_resp_potentialFP)

openxlsx::write.xlsx(vene_astral_resp_final, "Results/VenetoclaxAstral_RESP_final_summary.xlsx")

# proteins reported in Ramos et al.
c("ARHGAP25", "DDX46", "DKX1", "FHOD1", "IWS1", "NOP14", 
  "PRPF40A", "SF3B1", "ZMYND8")

sum(!(unique(vene_astral_resp_final$Gene[which(vene_astral_resp_final$kept_for_further_analysis & 
                                                 !vene_astral_resp_final$in_merCASBA)]) %in%
        c("ARHGAP25", "DDX46", "DKX1", "FHOD1", "IWS1", "NOP14", 
          "PRPF40A", "SF3B1", "ZMYND8"))
    ) # 52 never reported before
