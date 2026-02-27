library(IMPRINTS.CETSA.app)

veneastral_resp_final <- openxlsx::read.xlsx("Results/VenetoclaxAstral_RESP_final_summary.xlsx")

mercasba <- openxlsx::read.xlsx("Data/mapped_merCASBA.xlsx")

unique(veneastral_resp_final$Gene[veneastral_resp_final$Gene %in% unique(mercasba$Gene)])
# PARP1 BCAP31 LMNB1 (MATR3)
unique(veneastral_resp_final$Gene[!(veneastral_resp_final$Gene %in% unique(mercasba$Gene))])
# ZMYND8 ARHGAP25 IWS1 
# ATF7IP CTDP1 SF3B1        
                       
veneastral_pepnorm <- readr::read_tsv("Data/VenetoclaxAstral/VenetoclaxAstral_peptides_postnorm.txt")

### old
v1_3h <- imprints_sequence_peptides(veneastral_pepnorm[,-grep("12h$", colnames(veneastral_pepnorm))], control = "Vehicle3h",
                                    dataset_name = NULL,
                                    proteins = c("P09874", "P51572", "P20700", "P43243"),
                                    sequence = c("209-221", "160-167", "220-241", "638-701")) 
v1_3h <- imprints_remove_peptides(v1_3h,
                                  proteins = c("P09874", "P51572", "P20700", "P43243"),
                                  sequence = c("209-221", "160-167", "220-241", "638-701")) 
v1_3h <- v1_3h[,-grep("Vehicle", colnames(v1_3h))]

v1_12h <- imprints_sequence_peptides(veneastral_pepnorm[,-grep("3h$", colnames(veneastral_pepnorm))], control = "Vehicle12h",
                                    dataset_name = NULL,
                                    proteins = c("P09874", "P51572", "P20700", "P43243"),
                                    sequence = c("209-221", "160-167", "220-241", "638-701")) 
v1_12h <- imprints_remove_peptides(v1_12h,
                                  proteins = c("P09874", "P51572", "P20700", "P43243"),
                                  sequence = c("209-221", "160-167", "220-241", "638-701")) 
v1_12h <- v1_12h[,-grep("Vehicle", colnames(v1_12h))]

v1 <- imprints_join_peptides(list(v1_3h, v1_12h), "exact")


p <- imprints_barplotting_peptides(v1, treatmentlevel = c("Venetoclax3h", "Venetoclax12h"),
                              colorpanel = c("#ff8080", "#d090f8"),
                              format = "RESP_peptide")
ggpubr::ggarrange(plotlist = p[-3], common.legend = TRUE, legend = "right",
                  nrow = 1, ncol = 3)
ggsave("KnownCaspaseTarget_AstralVenetoclax.png", dpi = 300, width = 12, height = 8)



### new
v2_3h <- imprints_sequence_peptides(veneastral_pepnorm[,-grep("12h$", colnames(veneastral_pepnorm))], control = "Vehicle3h",
                                    dataset_name = NULL,
                                    proteins = c("Q9ULU4", "P42331", "Q96ST2", 
                                                 "O75533", "Q9Y5B0", "Q6VMQ6"),
                                    sequence = c("674-950", "356-377", "447-566", 
                                                 "271-300", "800-851", "341-369")) 
v2_3h <- imprints_remove_peptides(v2_3h,
                                  proteins = c("Q9ULU4", "P42331", "Q96ST2", 
                                               "O75533", "Q9Y5B0", "Q6VMQ6"),
                                  sequence = c("674-950", "356-377", "447-566", 
                                               "271-300", "800-851", "341-369")) 
v2_3h <- v2_3h[,-grep("Vehicle", colnames(v2_3h))]

v2_12h <- imprints_sequence_peptides(veneastral_pepnorm[,-grep("3h$", colnames(veneastral_pepnorm))], control = "Vehicle12h",
                                     dataset_name = NULL,
                                     proteins = c("Q9ULU4", "P42331", "Q96ST2", 
                                                  "O75533", "Q9Y5B0", "Q6VMQ6"),
                                     sequence = c("674-950", "356-377", "447-566", 
                                                  "271-300", "800-851", "341-369")) 
v2_12h <- imprints_remove_peptides(v2_12h,
                                   proteins = c("Q9ULU4", "P42331", "Q96ST2", 
                                                "O75533", "Q9Y5B0", "Q6VMQ6"),
                                   sequence = c("674-950", "356-377", "447-566", 
                                                "271-300", "800-851", "341-369")) 
v2_12h <- v2_12h[,-grep("Vehicle", colnames(v2_12h))]

v2 <- imprints_join_peptides(list(v2_3h, v2_12h), "exact")


p <- imprints_barplotting_peptides(v2, treatmentlevel = c("Venetoclax3h", "Venetoclax12h"),
                                   colorpanel = c("#ff8080", "#d090f8"),
                                   format = "RESP_peptide")

ggpubr::ggarrange(plotlist = p[c(2,4,5)], common.legend = TRUE, legend = "right",
                  nrow = 1, ncol = 3)
ggsave("Figures/Fig2C-Unknown1CaspaseTarget_AstralVenetoclax.png", dpi = 300, width = 12, height = 8)

ggpubr::ggarrange(plotlist = p[c(3,1,6)], common.legend = TRUE, legend = "right",
                  nrow = 1, ncol = 3)
ggsave("Figures/Fig2C-Unknown2CaspaseTarget_AstralVenetoclax.png", dpi = 300, width = 12, height = 8)


rm(veneastral_resp_final, mercasba, veneastral_pepnorm, 
   v1_3h, v1_12h, v1, v2_3h, v2_12h, v2, p)
