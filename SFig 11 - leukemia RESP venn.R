gemO19 <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_Gemcitabine_OCILY19",
                                         pattern = "full", full.names = TRUE))
gemO3 <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_Gemcitabine_OCILY3",
                                        pattern = "full", full.names = TRUE))
atiap4h <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_AT-IAP4h",
                                          pattern = "summary", full.names = TRUE))

astral <- openxlsx::read.xlsx("Results/VenetoclaxAstral_RESP_final_summary.xlsx")

# print venn diagram 
venn <- VennDiagram::venn.diagram(list("Kasumi-1\nVenetoclax" = unique(astral$Gene[which(astral$kept_for_further_analysis)]),
                                       "MOLM-16\nAT-IAP 4h" = atiap4h$Gene[which(atiap4h$category == "RESP")],
                                       "OCI-LY3\nGemcitabine 8h" = gemO3$Gene[which(gemO3$category == "RESP")],
                                       "OCI-LY19\nGemcitabine 8h" = gemO19$Gene[which(gemO19$category == "RESP")]
                                       ), 
                                  filename = NULL,
                                  disable.logging = TRUE, margin = 0.5,
                                  fill = c("purple", "#fca311", "#94d2bd", "#005f73"), lwd = 2, cex = 3,
                                  main = "Venn Diagram of proteins with\na significant RESP effect\nin different IMPRINTS-CETSA datasets",
                                  main.fontface = "bold", main.cex = 3, main.pos = c(0.5,0.8),
                                  #cat.dist = c(0.06,0.0675),
                                  fontface = "bold", fontfamiliy = "sans",
                                  cat.cex = 2.4, cat.fontface = "bold",
                                  cat.fontfamily = "sans")

png("Figures/RESP_leukemia_VennDiagram.png", 
    width=16, height=16, units = "in", res=300)
grid::grid.newpage()
grid::grid.draw(venn)
dev.off()
