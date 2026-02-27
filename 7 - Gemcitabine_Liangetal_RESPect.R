library(IMPRINTS.CETSA.app)


### Gemcitabine OCILY3 1h / 8h
gemoc3_norm <- readr::read_tsv("Data/Gemcitabine_Liangetal/OCILY3/GemcitabineOCILY3_peptides_post_norm.txt")

# compute log2 fold-changes
gemoc3_diff <- imprints_sequence_peptides(gemoc3_norm, control = "GemO318.OCILY3.Vehicle.1h", 
                                          dataset_name = "GemcitabineOCILY3_peptides_foldchange")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_CaldiffPeptides_GemcitabineOCILY3_peptides_foldchange"),
            "Results/GemcitabineOCILY3_peptides_foldchange.txt")

# run RESPect
gemoc3_respect <- imprints_cleaved_peptides(gemoc3_norm, gemoc3_diff, 
                                            control = "GemO318.OCILY3.Vehicle.1h", 
                                            folder_name = "Gemcitabine_OCILY3",
                                            fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                            FDR = 0.01, comb_pv_method = "george",
                                            min_peptide = 5, min_ValidValue = 0.5)

# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_Gemcitabine_OCILY3"),
            "Results/RESP_analysis_Gemcitabine_OCILY3")

# remake volcano plot from RESPect Gemcitabine OCILY3 - SFig 9
gemoc3_respect_full <- openxlsx::read.xlsx(list.files(path = "RESP_analysis_Gemcitabine_OCILY3",
                                                     pattern = "full", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "RESP_analysis_Gemcitabine_OCILY3",
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 2)
cutoff <- cutoff[which(cutoff$treatment == "GemO318.OCILY3.Gemcitabine.8h"),]

gemoc3_respect_full <- gemoc3_respect_full[which(gemoc3_respect_full$treatment == "GemO318.OCILY3.Gemcitabine.8h"),1:7]
gemoc3_respect_full$treatment <- "OCI-LY3\nGemcitabine 8h vs DMSO"
cutoff$treatment <- "OCI-LY3\nGemcitabine 8h vs DMSO"
gemoc3_respect_full <- left_join(gemoc3_respect_full, gemoc3_respect[,c("id", "Gene", "description", "category")], by = c("id", "Gene", "description"))
gemoc3_respect_full$category <- sub("m$", "", gemoc3_respect_full$category)
gemoc3_respect_full$category[which(is.na(gemoc3_respect_full$category))] <- "Non Significant"
gemoc3_respect <- mapply(function(t, n){
  paste0(n, " (n=", t, ")")
}, table(gemoc3_respect_full$category), names(table(gemoc3_respect_full$category)))
gemoc3_respect_full$category <- unname(gemoc3_respect[gemoc3_respect_full$category])

gemoc3_respect_full$category <- factor(gemoc3_respect_full$category, 
                                      levels = c(paste0("RESP (n=", sum(grepl("^RESP ", gemoc3_respect_full$category)), ")"),
                                                 paste0("SP (n=", sum(grepl("^SP ", gemoc3_respect_full$category)), ")"),
                                                 paste0("MP (n=", sum(grepl("^MP ", gemoc3_respect_full$category)), ")"),
                                                 paste0("Unconfident (n=", sum(grepl("^Unconfident ", gemoc3_respect_full$category)), ")"),
                                                 paste0("Non Significant (n=", sum(grepl("^Non ", gemoc3_respect_full$category)), ")"))
                                      )


df_curve <- data.frame(RESP_score = seq(-max(abs(gemoc3_respect_full$RESP_score), na.rm = TRUE),
                                        max(abs(gemoc3_respect_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "OCI-LY3\nGemcitabine 8h vs DMSO"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, -0.5, 0.5, cutoff$pval))

xadd <- as.data.frame(matrix(nrow = 3, ncol = ncol(gemoc3_respect_full), 
                             dimnames = list(c(1:3), colnames(gemoc3_respect_full))
                             )
                      ) 
xadd$treatment <- "OCI-LY3\nGemcitabine 8h vs DMSO"

if(length(grep("n=0", levels(gemoc3_respect_full$category), value = TRUE))){
  xadd$category <- grep("n=0", levels(gemoc3_respect_full$category), value = TRUE)
  gemoc3_respect_full <- as.data.frame(rbind(gemoc3_respect_full, xadd))
}

ggplot(gemoc3_respect_full, aes(RESP_score, -log10(combined_pvalue), color = category)) +
  geom_point() +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(gemoc3_respect_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(gemoc3_respect_full$RESP_score), na.rm = TRUE), max(abs(gemoc3_respect_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESP Volcano plot",
       y = "-log10(Combined p-value)", x = "RESP score") +
  scale_color_manual(values = rlang::set_names(c("blue", "#F7BD00", "#F76000", "red", "grey70"), 
                                               levels(gemoc3_respect_full$category))) +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = gemoc3_respect_full[grep("^RESP ", gemoc3_respect_full$category),], 
                            aes(RESP_score, -log10(combined_pvalue), label = Gene), 
                            show.legend = FALSE, min.segment.length = 0, size = 6) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", color = "black", size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))

ggsave("Figures/RESPVolcano_GemOCI-LY3_8h.png", dpi = 300, device = "png", width = 16, height = 12)


### Gemcitabine OCILY19 1h / 8h
gemoc19_norm <- readr::read_tsv("Data/Gemcitabine_Liangetal/OCILY19/GemcitabineOCILY19_peptides_post_norm.txt")

# compute log2 fold-changes
gemoc19_diff <- imprints_sequence_peptides(gemoc19_norm, control = "GemO1918.OCILY19.Vehicle.1h", 
                                           dataset_name = "GemcitabineOCILY19_peptides_foldchange")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_CaldiffPeptides_GemcitabineOCILY19_peptides_foldchange"),
            "Results/GemcitabineOCILY19_peptides_foldchange.txt")

gemoc19_respect <- imprints_cleaved_peptides(gemoc19_norm, gemoc19_diff, 
                                             control = "GemO1918.OCILY19.Vehicle.1h", 
                                             folder_name = "Gemcitabine_OCILY19",
                                             fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                             FDR = 0.01, comb_pv_method = "george",
                                             min_peptide = 5, min_ValidValue = 0.5)
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_Gemcitabine_OCILY19"),
            "Results/RESP_analysis_Gemcitabine_OCILY19")

# remake volcano plot from RESPect Gemcitabine OCILY19 - SFig 9
gemoc19_respect_full <- openxlsx::read.xlsx(list.files(path = "RESP_analysis_Gemcitabine_OCILY19",
                                                      pattern = "full", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "RESP_analysis_Gemcitabine_OCILY19",
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 2)
cutoff <- cutoff[which(cutoff$treatment == "GemO1918.OCILY19.Gemcitabine.8h"),]

gemoc19_respect_full <- gemoc19_respect_full[which(gemoc19_respect_full$treatment == "GemO1918.OCILY19.Gemcitabine.8h"),1:7]
gemoc19_respect_full$treatment <- "OCI-LY19\nGemcitabine 8h vs DMSO"
cutoff$treatment <- "OCI-LY19\nGemcitabine 8h vs DMSO"
gemoc19_respect_full <- left_join(gemoc19_respect_full, gemoc19_respect[,c("id", "Gene", "description", "category")], by = c("id", "Gene", "description"))
gemoc19_respect_full$category <- sub("m$", "", gemoc19_respect_full$category)
gemoc19_respect_full$category[which(is.na(gemoc19_respect_full$category))] <- "Non Significant"
gemoc19_respect <- mapply(function(t, n){
  paste0(n, " (n=", t, ")")
}, table(gemoc19_respect_full$category), names(table(gemoc19_respect_full$category)))
gemoc19_respect_full$category <- unname(gemoc19_respect[gemoc19_respect_full$category])

gemoc19_respect_full$category <- factor(gemoc19_respect_full$category, 
                                       levels = c(paste0("RESP (n=", sum(grepl("^RESP ", gemoc19_respect_full$category)), ")"),
                                                  paste0("SP (n=", sum(grepl("^SP ", gemoc19_respect_full$category)), ")"),
                                                  paste0("MP (n=", sum(grepl("^MP ", gemoc19_respect_full$category)), ")"),
                                                  paste0("Unconfident (n=", sum(grepl("^Unconfident ", gemoc19_respect_full$category)), ")"),
                                                  paste0("Non Significant (n=", sum(grepl("^Non ", gemoc19_respect_full$category)), ")"))
                                       )


df_curve <- data.frame(RESP_score = seq(-max(abs(gemoc19_respect_full$RESP_score), na.rm = TRUE),
                                        max(abs(gemoc19_respect_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "OCI-LY19\nGemcitabine 8h vs DMSO"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, -0.5, 0.5, cutoff$pval))

xadd <- as.data.frame(matrix(nrow = 3, ncol = ncol(gemoc19_respect_full), 
                             dimnames = list(c(1:3), colnames(gemoc19_respect_full))
                             )
                      ) 
xadd$treatment <- "OCI-LY19\nGemcitabine 8h vs DMSO"

if(length(grep("n=0", levels(gemoc19_respect_full$category), value = TRUE))){
  xadd$category <- grep("n=0", levels(gemoc19_respect_full$category), value = TRUE)
  gemoc19_respect_full <- as.data.frame(rbind(gemoc19_respect_full, xadd))
}

ggplot(gemoc19_respect_full, aes(RESP_score, -log10(combined_pvalue), color = category)) +
  geom_point() +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(gemoc19_respect_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(gemoc19_respect_full$RESP_score), na.rm = TRUE), max(abs(gemoc19_respect_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESP Volcano plot",
       y = "-log10(Combined p-value)", x = "RESP score") +
  scale_color_manual(values = rlang::set_names(c("blue", "#F7BD00", "#F76000", "red", "grey70"), 
                                               levels(gemoc19_respect_full$category))) +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = gemoc19_respect_full[grep("^RESP ", gemoc19_respect_full$category),], 
                            aes(RESP_score, -log10(combined_pvalue), label = Gene), 
                            show.legend = FALSE, min.segment.length = 0, size = 6) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", color = "black", size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))

ggsave("Figures/RESPVolcano_GemOCI-LY19_8h.png", dpi = 300, device = "png", width = 16, height = 12)

rm(gemoc3_respect_full, gemoc19_respect_full, cutoff, df_cuvre, xadd)
