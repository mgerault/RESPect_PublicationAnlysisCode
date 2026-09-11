library(IMPRINTS.CETSA.app)

### Replotting Volcano plot from hits from IMPRINTS-CETSA protein data
# data shaping venetoclax 
x1 <- openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_3h_hits_analysis",
                                     pattern = "analysis_tab", full.names = TRUE)
                          )
x1 <- x1[,c("id", "Gene", grep("Combinedpval|IS_|category", colnames(x1), value = TRUE))]
colnames(x1) <- sub("_.*", "", colnames(x1))
x1$treatment <- "Venetoclax 3h"

x2 <- openxlsx::read.xlsx(list.files(path = "Results/venetoclax_astral_12h_hits_analysis",
                                     pattern = "analysis_tab", full.names = TRUE))
x2 <- x2[,c("id", "Gene", grep("Combinedpval|IS_|category", colnames(x2), value = TRUE))]
colnames(x2) <- sub("_.*", "", colnames(x2))
x2$treatment <- "Venetoclax 12h"

x <- as.data.frame(rbind(x1, x2))

cutoff1 <- readr::read_tsv(list.files(path = "Results/venetoclax_astral_3h_hits_analysis",
                                      pattern = "cutoff", full.names = TRUE),
                           n_max = 1)
cutoff2 <- readr::read_tsv(list.files(path = "Results/venetoclax_astral_12h_hits_analysis",
                                      pattern = "cutoff", full.names = TRUE),
                           n_max = 1)
cutoff1$treatment <- "Venetoclax 3h"
cutoff2$treatment <- "Venetoclax 12h"
cutoff <- as.data.frame(rbind(cutoff1, cutoff2))


df_curve <- as.data.frame(rbind(data.frame(IS = seq(-max(abs(x$IS), na.rm = TRUE),
                                                    max(abs(x$IS), na.rm = TRUE), 0.001),
                                           treatment = "Venetoclax 3h"),
                                data.frame(IS = seq(-max(abs(x$IS), na.rm = TRUE),
                                                    max(abs(x$IS), na.rm = TRUE), 0.001),
                                           treatment = "Venetoclax 12h")
                                )
                          )
df_curve <- df_curve %>% group_by(treatment) %>%
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(IS, -1.5, 1.5, 
                                                   cutoff$pval[cutoff$treatment == treatment[1]],
                                                   curvature = 0.05)
                )

x$hits <- x$category != "NN"
x$treatment <- factor(x$treatment, levels = c("Venetoclax 3h", "Venetoclax 12h"))
df_curve$treatment <- factor(df_curve$treatment, levels = c("Venetoclax 3h", "Venetoclax 12h"))

ggplot(x, aes(IS, -log10(Combinedpval), color = hits)) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  geom_line(data = df_curve, aes(x = IS, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(x$Combinedpval), na.rm = TRUE))) +
  xlim(c(-max(abs(x$IS), na.rm = TRUE), max(abs(x$IS), na.rm = TRUE))) +
  labs(title = "",
       y = "-log10(Combined p-value)", x = "I-score") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = x[which(x$hits),], 
                            aes(IS, -log10(Combinedpval), label = Gene), 
                            show.legend = FALSE, min.segment.length = 0,
                            size = 6) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", color = "black", size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))

ggsave("Figures/ISVolcano_VenetoclaxAstral.png", dpi = 300, 
       width = 24, height = 14)

# data cleaning
rm(x, x1, x2, cutoff, cutoff1, cutoff2, df_curve)
