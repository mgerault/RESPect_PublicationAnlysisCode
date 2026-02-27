library(IMPRINTS.CETSA.app)

### AT-IAP 4h
# read IMPRINTS-CETSA post norm file
atiap_norm <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/AT-IAP4h_peptides_postnorm.txt")
# removing unnecessary treatment
atiap_norm <- atiap_norm[,-grep("_Vehicle", colnames(atiap_norm))]

# compute log2 fold-changes
atiap_diff <- imprints_sequence_peptides(atiap_norm, control = "TNFa4h", 
                                         dataset_name = "AT-IAP4h_peptides_foldchange")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_CaldiffPeptides_AT-IAP4h_peptides_foldchange"),
            "Results/AT-IAP4h_peptides_foldchange.txt")

# run RESPect
atiap_respect <- imprints_cleaved_peptides(atiap_norm, atiap_diff, 
                                           control = "TNFa4h", folder_name = "AT-IAP4h",
                                           fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                           FDR = 0.01, comb_pv_method = "george",
                                           min_peptide = 5, min_ValidValue = 0.5)
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_AT-IAP4h"),
            "Results/RESP_analysis_AT-IAP4h")


# remake volcano plot from RESPect AT-IAP 4h - SFig 9
atiap_respect_full <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_AT-IAP4h",
                                                     pattern = "full", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "Results/RESP_analysis_AT-IAP4h",
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 1)

atiap_respect_full <- atiap_respect_full[,1:7]
atiap_respect_full$treatment <- "MOLM-16\nAT-IPA 4h vs TNFa4h"
cutoff$treatment <- "MOLM-16\nAT-IPA 4h vs TNFa4h"
atiap_respect_full <- left_join(atiap_respect_full, atiap_respect[,c("id", "Gene", "description", "category")], by = c("id", "Gene", "description"))
atiap_respect_full$category <- sub("m$", "", atiap_respect_full$category)
atiap_respect_full$category[which(is.na(atiap_respect_full$category))] <- "Non Significant"
atiap_respect <- mapply(function(t, n){
  paste0(n, " (n=", t, ")")
}, table(atiap_respect_full$category), names(table(atiap_respect_full$category)))
atiap_respect_full$category <- unname(atiap_respect[atiap_respect_full$category])

atiap_respect_full$category <- factor(atiap_respect_full$category, 
                                      levels = c(paste0("RESP (n=", sum(grepl("^RESP ", atiap_respect_full$category)), ")"),
                                                 paste0("SP (n=", sum(grepl("^SP ", atiap_respect_full$category)), ")"),
                                                 paste0("MP (n=", sum(grepl("^MP ", atiap_respect_full$category)), ")"),
                                                 paste0("Unconfident (n=", sum(grepl("^Unconfident ", atiap_respect_full$category)), ")"),
                                                 paste0("Non Significant (n=", sum(grepl("^Non ", atiap_respect_full$category)), ")"))
                                      )


df_curve <- data.frame(RESP_score = seq(-max(abs(atiap_respect_full$RESP_score), na.rm = TRUE),
                                        max(abs(atiap_respect_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "MOLM-16\nAT-IPA 4h vs TNFa4h"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, -0.5, 0.5, cutoff$pval))

xadd <- as.data.frame(matrix(nrow = 3, ncol = ncol(atiap_respect_full), 
                             dimnames = list(c(1:3), colnames(atiap_respect_full))
                             )
                      ) 
xadd$treatment <- "MOLM-16\nAT-IPA 4h vs TNFa4h"

if(length(grep("n=0", levels(atiap_respect_full$category), value = TRUE))){
  xadd$category <- grep("n=0", levels(atiap_respect_full$category), value = TRUE)
  atiap_respect_full <- as.data.frame(rbind(atiap_respect_full, xadd))
}

ggplot(atiap_respect_full, aes(RESP_score, -log10(combined_pvalue), color = category)) +
  geom_point() +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(atiap_respect_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(atiap_respect_full$RESP_score), na.rm = TRUE), max(abs(atiap_respect_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESP Volcano plot",
       y = "-log10(Combined p-value)", x = "RESP score") +
  scale_color_manual(values = rlang::set_names(c("blue", "#F7BD00", "#F76000", "red", "grey70"), 
                                               levels(atiap_respect_full$category))) +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = atiap_respect_full[c(grep("^RESP ", atiap_respect_full$category),
                                                        which(!is.na(match(atiap_respect_full$Gene, c("DDX21")))
                                                              )
                                                        ),], 
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

ggsave("Figures/RESPVolcano_AT-IAP4h.png", dpi = 300, device = "png", width = 16, height = 12)

rm(atiap_respect_full, cutoff, df_cuvre, xadd)
