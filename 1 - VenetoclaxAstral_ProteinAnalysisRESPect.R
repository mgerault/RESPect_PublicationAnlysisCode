library(IMPRINTS.CETSA.app)
library(IMPRINTS.CETSA)

#### Protein level
vene_astral_protnorm <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_proteins_postnorm.txt")

vene_astral_protnorm3h <- vene_astral_protnorm[,-grep("12h", colnames(vene_astral_protnorm))]
vene_astral_protnorm12h <- vene_astral_protnorm[,-grep("3h", colnames(vene_astral_protnorm))]

vene_astral_protdiff3h <- imprints_caldiff_f(vene_astral_protnorm3h, reftreatment = "Vehicle3h")
# move file to Results folder
file.rename(list.files(path = list.files(pattern = "vene_astral_protnorm3h_\\d{6}_\\d{4}"), 
                       full.names = TRUE),
            "Results/VenetoclaxAstral_proteins_foldchange3h.txt")
unlink(list.files(pattern = "vene_astral_protnorm3h_\\d{6}_\\d{4}"), recursive = TRUE)

vene_astral_protdiff12h <- imprints_caldiff_f(vene_astral_protnorm12h, reftreatment = "Vehicle12h")
# move file to Results folder
file.rename(list.files(path = list.files(pattern = "vene_astral_protnorm12h_\\d{6}_\\d{4}"), 
                       full.names = TRUE),
            "Results/VenetoclaxAstral_proteins_foldchange12h.txt"
            )
unlink(list.files(pattern = "vene_astral_protnorm12h_\\d{6}_\\d{4}"), recursive = TRUE)

## hits analysis
vene_astral_prothits3h <- imprints_IS(vene_astral_protnorm3h, vene_astral_protdiff3h, 
                                      ctrl = "Vehicle3h", 
                                      IS_cutoff = 1.5, FDR = 0.01,
                                      fixed_score_cutoff = TRUE, 
                                      folder_name = "venetoclax_astral_3h_hits_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_3h_hits_analysis"),
            "Results/venetoclax_astral_3h_hits_analysis")

vene_astral_prothits12h <- imprints_IS(vene_astral_protnorm12h, vene_astral_protdiff12h, 
                                       ctrl = "Vehicle12h", 
                                       IS_cutoff = 1.5, FDR = 0.01,
                                       fixed_score_cutoff = TRUE,
                                       folder_name = "venetoclax_astral_12h_hits_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_12h_hits_analysis"),
            "Results/venetoclax_astral_12h_hits_analysis")

### ORA on CETSA hits (SFig 5)
vene_astral_protORA_GOBP <- compare_enrich(as.data.frame(rbind(vene_astral_prothits3h, vene_astral_prothits12h)) %>%
                 mutate(treatment = factor(treatment, 
                                           levels = c("Venetoclax3h", "Venetoclax12h"))
                        ), treatment_column = "treatment",
               pval_cutoff = 0.01, minGSSize = 3,
               database = "GO", ont = "BP")

ggsave("Figures/CompareEnrichPlot_GO-BP_AstralVenetoclax_IMPRINTShits.png",
       plot = vene_astral_protORA_GOBP$graph,
       dpi = 300, width = 8, height = 8)

vene_astral_protORA_GOCC <- compare_enrich(as.data.frame(rbind(vene_astral_prothits3h, vene_astral_prothits12h)) %>%
                 mutate(treatment = factor(treatment, 
                                           levels = c("Venetoclax3h", "Venetoclax12h"))
                        ), treatment_column = "treatment",
               pval_cutoff = 0.01, minGSSize = 3,
               database = "GO", ont = "CC")

ggsave("Figures/CompareEnrichPlot_GO-CC_AstralVenetoclax_IMPRINTShits.png",
       plot = vene_astral_protORA_GOCC$graph,
       dpi = 300, width = 8, height = 8)

## hits analysis QP
# Stringent: FDR 5% + |FC| >= 0.5
vene_astral_protqphits3h <- IMPRINTS.PhosphoQP::imprints_phoQP_hit_protein(vene_astral_protnorm3h, 
                                                                           vene_astral_protdiff3h, 
                                                                           ctrl = "Vehicle3h", QP_name = "36C",
                                                                           FC_cutoff = 0.5, FDR = 0.05, curvature = 0.05,
                                                                           folder_name = "venetoclax_astral_3h_QPHits_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_3h_QPHits_analysis"),
            "Results/venetoclax_astral_3h_QPHits_analysis")

vene_astral_protqphits12h <- IMPRINTS.PhosphoQP::imprints_phoQP_hit_protein(vene_astral_protnorm12h, 
                                                                            vene_astral_protdiff12h, 
                                                                            ctrl = "Vehicle12h", QP_name = "36C",
                                                                            FC_cutoff = 0.5, FDR = 0.05, curvature = 0.05,
                                                                            folder_name = "venetoclax_astral_12h_QPHits_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_12h_QPHits_analysis"),
            "Results/venetoclax_astral_12h_QPHits_analysis")


## analyzing potential DEAR effect: QP vs 37°C
# 3 hours
vene_astral_DEAR3h <- vene_astral_protdiff3h[,-grep("^5", colnames(vene_astral_protdiff3h))]
colnames(vene_astral_DEAR3h)[grep("^36C_", colnames(vene_astral_DEAR3h))] <- paste0("11C",
                                                                                    sub("36C", "", 
                                                                                        sub("Venetoclax3h", "QP-Venetoclax3h", colnames(vene_astral_DEAR3h)[grep("^36C_", colnames(vene_astral_DEAR3h))])
                                                                                        )
                                                                                    )
colnames(vene_astral_DEAR3h)[grep("^37C_", colnames(vene_astral_DEAR3h))] <- paste0("11C",
                                                                                    sub("37C", "", 
                                                                                        sub("Venetoclax3h", "37C-Venetoclax3h", colnames(vene_astral_DEAR3h)[grep("^37C_", colnames(vene_astral_DEAR3h))])
                                                                                        )
                                                                                    )
# computing QP - 37°C
vene_astral_DEAR3hdiff3h <- imprints_caldiff_f(vene_astral_DEAR3h, "37C-Venetoclax3h")
# statistical analysis
vene_astral_DEAR_hits3h <- IMPRINTS.PhosphoQP::imprints_phoQP_hit_protein(vene_astral_DEAR3h, 
                                                                          vene_astral_DEAR3hdiff3h,
                                                                          ctrl = "37C-Venetoclax3h", QP_name = "11C",
                                                                          FC_cutoff = 0.2, FDR = 0.05, curvature = 0.05,
                                                                          folder_name = "venetoclax_astral_3h_DEAR_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_3h_DEAR_analysis"),
            "Results/venetoclax_astral_3h_DEAR_analysis")

# 12 hours
vene_astral_DEAR12h <- vene_astral_protdiff12h[,-grep("^5", colnames(vene_astral_protdiff12h))]
colnames(vene_astral_DEAR12h)[grep("^36C_", colnames(vene_astral_DEAR12h))] <- paste0("11C",
                                                                                      sub("36C", "", 
                                                                                          sub("Venetoclax12h", "QP-Venetoclax12h", colnames(vene_astral_DEAR12h)[grep("^36C_", colnames(vene_astral_DEAR12h))])
                                                                                          )
                                                                                      )
colnames(vene_astral_DEAR12h)[grep("^37C_", colnames(vene_astral_DEAR12h))] <- paste0("11C",
                                                                                      sub("37C", "", 
                                                                                          sub("Venetoclax12h", "37C-Venetoclax12h", colnames(vene_astral_DEAR12h)[grep("^37C_", colnames(vene_astral_DEAR12h))])
                                                                                          )
                                                                                      )
# computing QP - 37°C
vene_astral_DEAR3hdiff12h <- imprints_caldiff_f(vene_astral_DEAR12h, "37C-Venetoclax12h")
# statistical analysis
vene_astral_DEAR_hits12h <- IMPRINTS.PhosphoQP::imprints_phoQP_hit_protein(vene_astral_DEAR12h, 
                                                                           vene_astral_DEAR3hdiff12h,
                                                                           ctrl = "37C-Venetoclax12h", QP_name = "11C",
                                                                           FC_cutoff = 0.2, FDR = 0.05, curvature = 0.05,
                                                                           folder_name = "venetoclax_astral_12h_DEAR_analysis")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_venetoclax_astral_12h_DEAR_analysis"),
            "Results/venetoclax_astral_12h_DEAR_analysis")


### Comparing number of quantified proteins with Astral MS compare to previous venetoclax dataset
pr1h <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax1h_proteins_postnorm.txt")
pr6h <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax6h_proteins_postnorm.txt")
pras <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_proteins_postnorm.txt")

pep1h <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax1h_peptides_postnorm.txt")
pep6h <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax6h_peptides_postnorm.txt")
pepas <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_peptides_postnorm.txt")

nrow(pras) / c(nrow(pr1h), nrow(pr6h)) # 1.992623 1.935510
nrow(pepas) / c(nrow(pep1h), nrow(pep6h)) #  1.822965 1.950267
# almost double of proteins and peptides correctly quantified

rm(pr1h, pr6h, pras, pep1h, pep6h, pepas)


#### Peptide level
vene_astral_pepnorm <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_peptides_postnorm.txt")

vene_astral_pepnorm3h <- vene_astral_pepnorm[,-grep("12h", colnames(vene_astral_pepnorm))]
vene_astral_pepnorm12h <- vene_astral_pepnorm[,-grep("3h", colnames(vene_astral_pepnorm))]

vene_astral_pepdiff3h <- imprints_sequence_peptides(vene_astral_pepnorm3h, control = "Vehicle3h",
                                                    dataset_name = "VenetoclaxAstral3h")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_CaldiffPeptides_VenetoclaxAstral3h"),
            "Results/VenetoclaxAstral_peptides_foldchange_3h.txt")

vene_astral_pepdiff12h <- imprints_sequence_peptides(vene_astral_pepnorm12h, control = "Vehicle12h",
                                                     dataset_name = "VenetoclaxAstral12h")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_CaldiffPeptides_VenetoclaxAstral12h"),
            "Results/VenetoclaxAstral_peptides_foldchange_12h.txt")

## RESPect analysis
# 3h
vene_astral_resp3h <- imprints_cleaved_peptides(vene_astral_pepnorm3h, vene_astral_pepdiff3h, 
                                                control = "Vehicle3h", min_peptide = 5,
                                                min_ValidValue = 0.5, FDR = 0.005, RESP_score = 0.5,
                                                fixed_score_cutoff = TRUE, 
                                                folder_name = "venetoclax_astral_3h")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_venetoclax_astral_3h"),
            "Results/RESP_analysis_venetoclax_astral_3h")

table(vene_astral_resp3h$category)
# MP   MPm  RESP  SP  SPm Unconfident
# 10     3   70   17   3          21

# remake volcano plot from RESPect Venetoclax 3h - Fig 2
vene3h_respect <- vene_astral_resp3h
vene3h_respect_full <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_venetoclax_astral_3h",
                                                      pattern = "full", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "Results/RESP_analysis_venetoclax_astral_3h",
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 1)

vene3h_respect_full <- vene3h_respect_full[,1:7]
vene3h_respect_full$treatment <- "Venetoclax 3h"
cutoff$treatment <- "Venetoclax 3h"
vene3h_respect_full <- left_join(vene3h_respect_full, vene3h_respect[,c("id", "Gene", "description", "category")], by = c("id", "Gene", "description"))
vene3h_respect_full$category <- sub("m$", "", vene3h_respect_full$category)
vene3h_respect_full$category[which(is.na(vene3h_respect_full$category))] <- "Non Significant"
vene3h_respect <- mapply(function(t, n){
  paste0(n, " (n=", t, ")")
}, table(vene3h_respect_full$category), names(table(vene3h_respect_full$category)))
vene3h_respect_full$category <- unname(vene3h_respect[vene3h_respect_full$category])

vene3h_respect_full$category <- factor(vene3h_respect_full$category, 
                                      levels = c(paste0("RESP (n=", sum(grepl("^RESP ", vene3h_respect_full$category)), ")"),
                                                 paste0("SP (n=", sum(grepl("^SP ", vene3h_respect_full$category)), ")"),
                                                 paste0("MP (n=", sum(grepl("^MP ", vene3h_respect_full$category)), ")"),
                                                 paste0("Unconfident (n=", sum(grepl("^Unconfident ", vene3h_respect_full$category)), ")"),
                                                 paste0("Non Significant (n=", sum(grepl("^Non ", vene3h_respect_full$category)), ")"))
                                      )


df_curve <- data.frame(RESP_score = seq(-max(abs(vene3h_respect_full$RESP_score), na.rm = TRUE),
                                        max(abs(vene3h_respect_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "Venetoclax 3h"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, -0.5, 0.5, cutoff$pval))

xadd <- as.data.frame(matrix(nrow = 3, ncol = ncol(vene3h_respect_full), 
                             dimnames = list(c(1:3), colnames(vene3h_respect_full))
                             )
                      ) 
xadd$treatment <- "Venetoclax 3h"

if(length(grep("n=0", levels(vene3h_respect_full$category), value = TRUE))){
  xadd$category <- grep("n=0", levels(vene3h_respect_full$category), value = TRUE)
  vene3h_respect_full <- as.data.frame(rbind(vene3h_respect_full, xadd))
}

ggplot(vene3h_respect_full, aes(RESP_score, -log10(combined_pvalue), color = category)) +
  geom_point() +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(vene3h_respect_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(vene3h_respect_full$RESP_score), na.rm = TRUE), max(abs(vene3h_respect_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESP Volcano plot",
       y = "-log10(Combined p-value)", x = "RESP score") +
  scale_color_manual(values = rlang::set_names(c("blue", "#F7BD00", "#F76000", "red", "grey70"), 
                                               levels(vene3h_respect_full$category))) +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = vene3h_respect_full[grep("^RESP ", vene3h_respect_full$category),], 
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

ggsave("Figures/RESPVolcano_Venetoclax3h.png", dpi = 300, device = "png", width = 16, height = 12)

rm(vene3h_respect, vene3h_respect_full, cutoff, df_curve, xadd)

# Plot hits if necessary
#imprints_barplotting_categorize_peptides(vene_astral_pepnorm3h, vene_astral_resp3h,
#                                         "Venetoclax3h", "Vehicle3h", 
#                                         "peptide_one", "red",
#                                         "venetoclax_astral_3h_RESPect_barplotting")

# 12h
vene_astral_resp12h <- imprints_cleaved_peptides(vene_astral_pepnorm12h, vene_astral_pepdiff12h, 
                                                 control = "Vehicle12h", min_peptide = 5,
                                                 min_ValidValue = 0.5, 
                                                 FDR = 0.005, RESP_score = 0.5,
                                                 fixed_score_cutoff = TRUE, 
                                                 folder_name = "venetoclax_astral_12h")

# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_venetoclax_astral_12h"),
            "Results/RESP_analysis_venetoclax_astral_12h")

table(vene_astral_resp12h$category)
# MP   MPm  RESP  SP  SPm Unconfident
# 19     6   49   43    7         29


# remake volcano plot from RESPect Venetoclax 12h - Fig 2
vene12h_respect <- vene_astral_resp12h
vene12h_respect_full <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_venetoclax_astral_12h",
                                                       pattern = "full", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "Results/RESP_analysis_venetoclax_astral_12h",
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 1)

vene12h_respect_full <- vene12h_respect_full[,1:7]
vene12h_respect_full$treatment <- "Venetoclax 12h"
cutoff$treatment <- "Venetoclax 12h"
vene12h_respect_full <- left_join(vene12h_respect_full, vene12h_respect[,c("id", "Gene", "description", "category")], by = c("id", "Gene", "description"))
vene12h_respect_full$category <- sub("m$", "", vene12h_respect_full$category)
vene12h_respect_full$category[which(is.na(vene12h_respect_full$category))] <- "Non Significant"
vene12h_respect <- mapply(function(t, n){
  paste0(n, " (n=", t, ")")
}, table(vene12h_respect_full$category), names(table(vene12h_respect_full$category)))
vene12h_respect_full$category <- unname(vene12h_respect[vene12h_respect_full$category])

vene12h_respect_full$category <- factor(vene12h_respect_full$category, 
                                       levels = c(paste0("RESP (n=", sum(grepl("^RESP ", vene12h_respect_full$category)), ")"),
                                                  paste0("SP (n=", sum(grepl("^SP ", vene12h_respect_full$category)), ")"),
                                                  paste0("MP (n=", sum(grepl("^MP ", vene12h_respect_full$category)), ")"),
                                                  paste0("Unconfident (n=", sum(grepl("^Unconfident ", vene12h_respect_full$category)), ")"),
                                                  paste0("Non Significant (n=", sum(grepl("^Non ", vene12h_respect_full$category)), ")"))
                                       )


df_curve <- data.frame(RESP_score = seq(-max(abs(vene12h_respect_full$RESP_score), na.rm = TRUE),
                                        max(abs(vene12h_respect_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "Venetoclax 12h"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, -0.5, 0.5, cutoff$pval))

xadd <- as.data.frame(matrix(nrow = 3, ncol = ncol(vene12h_respect_full), 
                             dimnames = list(c(1:3), colnames(vene12h_respect_full))
                             )
                      ) 
xadd$treatment <- "Venetoclax 12h"

if(length(grep("n=0", levels(vene12h_respect_full$category), value = TRUE))){
  xadd$category <- grep("n=0", levels(vene12h_respect_full$category), value = TRUE)
  vene12h_respect_full <- as.data.frame(rbind(vene12h_respect_full, xadd))
}

ggplot(vene12h_respect_full, aes(RESP_score, -log10(combined_pvalue), color = category)) +
  geom_point() +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(vene12h_respect_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(vene12h_respect_full$RESP_score), na.rm = TRUE), max(abs(vene12h_respect_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESP Volcano plot",
       y = "-log10(Combined p-value)", x = "RESP score") +
  scale_color_manual(values = rlang::set_names(c("blue", "#F7BD00", "#F76000", "red", "grey70"), 
                                               levels(vene12h_respect_full$category))) +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = vene12h_respect_full[grep("^RESP ", vene12h_respect_full$category),], 
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

ggsave("Figures/RESPVolcano_Venetoclax12h.png", dpi = 300, device = "png", width = 16, height = 12)

rm(vene12h_respect, vene12h_respect_full, cutoff, df_curve, xadd)

# Plot if necessary
#imprints_barplotting_categorize_peptides(vene_astral_pepnorm12h, vene_astral_resp12h,
#                                         "Venetoclax12h", "Vehicle12h", 
#                                         "peptide_one", "purple",
#                                         "venetoclax_astral_12h_RESPect_barplotting")


## Mapping to alternative splicing forms
# 3h
vene_astral_isoform3h <- imprints_isoform_peptides(vene_astral_pepnorm3h, vene_astral_resp3h, "Vehicle3h",
                                                   "Data/20240619_UP000005640_9606.fasta",
                                                   xlsxname = "VenetoclaxAstral_RESP_3h_isoform_ambiguity")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_3h_isoform_ambiguity.xlsx"),
            "Results/VenetoclaxAstral_RESP_3h_isoform_ambiguity.xlsx")

imprints_plotting_isoform_peptides(vene_astral_pepnorm3h, vene_astral_isoform3h, 
                                   "Vehicle3h", "Venetoclax3h",
                                   0.5, ret_plot = FALSE, save_pdf = TRUE,
                                   pdfname = "VenetoclaxAstral_RESP_3h_isoform_ambiguity_plots")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_3h_isoform_ambiguity_plots.pdf"),
            "Results/VenetoclaxAstral_RESP_3h_isoform_ambiguity_plots.pdf")

# 12h
vene_astral_isoform12h <- imprints_isoform_peptides(vene_astral_pepnorm12h, vene_astral_resp12h, "Vehicle12h",
                                                       "Data/20240619_UP000005640_9606.fasta",
                                                       xlsxname = "VenetoclaxAstral_RESP_12h_isoform_ambiguity")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_12h_isoform_ambiguity.xlsx"),
            "Results/VenetoclaxAstral_RESP_12h_isoform_ambiguity.xlsx")

imprints_plotting_isoform_peptides(vene_astral_pepnorm12h, vene_astral_isoform12h, 
                                   "Vehicle12h", "Venetoclax12h",
                                   0.5, ret_plot = FALSE, save_pdf = TRUE,
                                   pdfname = "VenetoclaxAstral_RESP_12h_isoform_ambiguity_plots")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_12h_isoform_ambiguity_plots.pdf"),
            "Results/VenetoclaxAstral_RESP_12h_isoform_ambiguity_plots.pdf")

## Mapping to PTMs
# 3h
vene_astral_ptm3h <- imprints_ptms_peptides(vene_astral_pepnorm3h, vene_astral_resp3h, 
                                            "Vehicle3h", xlsxname = "VenetoclaxAstral_RESP_3h_PTMs_mapping")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_3h_PTMs_mapping.xlsx"),
            "Results/VenetoclaxAstral_RESP_3h_PTMs_mapping.xlsx")

# 12h
vene_astral_ptm12h <- imprints_ptms_peptides(vene_astral_pepnorm12h, vene_astral_resp12h, 
                                             "Vehicle12h", xlsxname = "VenetoclaxAstral_RESP_12h_PTMs_mapping")
# move file to Results folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_VenetoclaxAstral_RESP_12h_PTMs_mapping.xlsx"),
            "Results/VenetoclaxAstral_RESP_12h_PTMs_mapping.xlsx")
