library(dplyr)
library(ggplot2)
library(ggraph)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(GOSemSim)

#### Volcanos miRNA
### Venetoclax 3H vs DMSO 3h --> Top 20
v3h_dmso3h_diff_mir <- readr::read_tsv("Data/RNAseq/sRNA/V3HvsDMSO3H.Differential_analysis_results.xls")
v3h_dmso3h_diff_mir$regulation <- ifelse(v3h_dmso3h_diff_mir$pval <= 0.05,
                                         ifelse(v3h_dmso3h_diff_mir$log2FoldChange > 1, "up",
                                                ifelse(v3h_dmso3h_diff_mir$log2FoldChange < -1, "down", "no")),
                                         "no")
v3h_dmso3h_diff_mir$regulation[which(is.na(v3h_dmso3h_diff_mir$regulation))] <- "no"

g <- ggplot(v3h_dmso3h_diff_mir, aes(log2FoldChange, -log10(pval),
                                 color = factor(regulation, levels = c("up", "down", "no"))
                                 )) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggrepel::geom_label_repel(data = v3h_dmso3h_diff_mir %>%
                              dplyr::filter(regulation == "up") %>%
                              mutate(rnk = order(order(-log10(pval)*log2FoldChange,
                                                       decreasing = TRUE))) %>%
                              dplyr::filter(rnk <= 20),
                            aes(log2FoldChange, -log10(pval), label = sRNA,
                                color = factor(regulation, levels = c("up", "down", "no"))),
                            min.segment.length = 0, show.legend = FALSE, size = 7.5, max.overlaps = 50) +
  ylim(c(0,4)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-6,6),
                     labels = paste(round(c(-6, -4, -2, 0, 2, 4, 6), 1))) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "no" = "grey30"),
                     labels = c("up" = paste0("up-regulated (n=", sum(v3h_dmso3h_diff_mir$regulation == "up"), ")"),
                                "down" = paste0("down-regulated (n=", sum(v3h_dmso3h_diff_mir$regulation == "down"), ")"),
                                "no" = paste0("non-significant (n=", sum(v3h_dmso3h_diff_mir$regulation == "no"), ")")),
                     name = "p-value <= 0.05\n|log2(fold-change)| >= 1") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold-change", y = "-log10(p-value)",
       title = "Venetoclax 3H vs DMSO 3H") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(size = 26),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(size = 26))

ggsave("volcano_V3HvsDMSO3H_miR_Top20.png", plot = g,
       device = "png", dpi = 300, height = 12, width = 16)


### Venetoclax 12h vs DMSO 12h --> Top 10
v12h_dmso12h_diff_mir <- readr::read_tsv("Data/RNAseq/sRNA/V12HvsDMSO12H.Differential_analysis_results.xls")
v12h_dmso12h_diff_mir$regulation <- ifelse(v12h_dmso12h_diff_mir$pval <= 0.05,
                                           ifelse(v12h_dmso12h_diff_mir$log2FoldChange > 1, "up",
                                                  ifelse(v12h_dmso12h_diff_mir$log2FoldChange < -1, "down", "no")),
                                           "no")
v12h_dmso12h_diff_mir$regulation[which(is.na(v12h_dmso12h_diff_mir$regulation))] <- "no"

g <- ggplot(v12h_dmso12h_diff_mir, aes(log2FoldChange, -log10(pval),
                                   color = factor(regulation, levels = c("up", "down", "no"))
                                   )) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggrepel::geom_label_repel(data = v12h_dmso12h_diff_mir %>%
                              dplyr::filter(regulation == "up") %>%
                              mutate(rnk = order(order(-log10(pval)*log2FoldChange,
                                                       decreasing = TRUE))) %>%
                              dplyr::filter(rnk <= 10),
                            aes(log2FoldChange, -log10(pval), label = sRNA,
                                color = factor(regulation, levels = c("up", "down", "no"))),
                            min.segment.length = 0, show.legend = FALSE, size = 7.5, max.overlaps = 50) +
  scale_x_continuous(breaks = c(-9, -6, -4, -2, 0, 2, 4, 6, 9), limits = c(-9,9),
                     labels = paste(round(c(-9, -6, -4, -2, 0, 2, 4, 6, 9), 1))) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "no" = "grey30"),
                     labels = c("up" = paste0("up-regulated (n=", sum(v12h_dmso12h_diff_mir$regulation == "up"), ")"),
                                "down" = paste0("down-regulated (n=", sum(v12h_dmso12h_diff_mir$regulation == "down"), ")"),
                                "no" = paste0("non-significant (n=", sum(v12h_dmso12h_diff_mir$regulation == "no"), ")")),
                     name = "p-value <= 0.05\n|log2(fold-change)| >= 1") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold-change", y = "-log10(p-value)",
       title = "Venetoclax 12H vs DMSO 12H") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(size = 26),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(size = 26))

ggsave("volcano_V12HvsDMSO12H_miR_Top10.png", plot = g,
       device = "png", dpi = 300, height = 12, width = 16)

rm(g)

#### Get most influential miRNA using Binomial Probability distribution and plot corresponding mRNA target
## Read data from mRNA level
v3h_dmso3h_diff_mrna <- openxlsx::read.xlsx("Data/RNAseq/mRNA/V3HvsDMSO3H.xlsx")
v3h_dmso3h_diff_mrna$regulation <- ifelse(v3h_dmso3h_diff_mrna$pvalue <= 0.05,
                                          ifelse(v3h_dmso3h_diff_mrna$log2FoldChange > 1, "up",
                                                 ifelse(v3h_dmso3h_diff_mrna$log2FoldChange < -1, "down", "no")),
                                          "no")
v3h_dmso3h_diff_mrna$regulation[which(is.na(v3h_dmso3h_diff_mrna$regulation))] <- "no"
colnames(v3h_dmso3h_diff_mrna)[1:2] <- c("target_gene", "target_gene_name")
v3h_dmso3h_diff_mrna$padj <- NULL

v12h_dmso12h_diff_mrna <- openxlsx::read.xlsx("Data/RNAseq/mRNA/V12HvsDMSO12H.xlsx")
v12h_dmso12h_diff_mrna$regulation <- ifelse(v12h_dmso12h_diff_mrna$pvalue <= 0.05,
                                            ifelse(v12h_dmso12h_diff_mrna$log2FoldChange > 1, "up",
                                                   ifelse(v12h_dmso12h_diff_mrna$log2FoldChange < -1, "down", "no")),
                                            "no")
v12h_dmso12h_diff_mrna$regulation[which(is.na(v12h_dmso12h_diff_mrna$regulation))] <- "no"
colnames(v12h_dmso12h_diff_mrna)[1:2] <- c("target_gene", "target_gene_name")
v12h_dmso12h_diff_mrna$padj <- NULL

## Read predicted gene target from miRNA
mir_predtarget_gene <- readr::read_tsv("Data/RNAseq/sRNA/targets_mRNA_gene.xls")
target_gene_name <- bitr(
  mir_predtarget_gene$target_gene,
  fromType = "ENSEMBL",
  toType = c("SYMBOL"),
  OrgDb = org.Hs.eg.db,
  drop = FALSE)
colnames(target_gene_name) <- c("target_gene", "target_gene_name")
mir_predtarget_gene <- left_join(mir_predtarget_gene, target_gene_name, by = "target_gene", relationship = "many-to-many")
colnames(mir_predtarget_gene)[1] <- "sRNA"

## summarise up-regulated miRNA
mir_up <- v3h_dmso3h_diff_mir[which(v3h_dmso3h_diff_mir$regulation == "up"),] %>%
  dplyr::select(sRNA, log2FoldChange, pval) %>%
  mutate(treatment = "Venetoclax 3H vs DMSO 3H") %>%
  rbind(v12h_dmso12h_diff_mir[which(v12h_dmso12h_diff_mir$regulation == "up"),] %>%
          dplyr::select(sRNA, log2FoldChange, pval) %>%
          mutate(treatment = "Venetoclax 12H vs DMSO 12H")
        ) %>%
  left_join(unique(mir_predtarget_gene[,c("sRNA", "target_gene", "target_gene_name")]), by = "sRNA",
            relationship = "many-to-many")

mir_up_mapped_mrna <- mir_up %>%
  dplyr::rename(treatment_mir = treatment, log2FoldChange_mir = log2FoldChange,
                pvalue_mir = pval) %>%
  left_join(v3h_dmso3h_diff_mrna %>%
              mutate(treatment_mrna = "Venetoclax 3H vs DMSO 3H") %>%
              dplyr::rename(log2FoldChange_mrna = log2FoldChange, pvalue_mrna = pvalue,
                            regulation_mrna = regulation) %>%
              rbind(v12h_dmso12h_diff_mrna %>%
                      mutate(treatment_mrna = "Venetoclax 12H vs DMSO 12H") %>%
                      dplyr::rename(log2FoldChange_mrna = log2FoldChange, pvalue_mrna = pvalue,
                                    regulation_mrna = regulation)
                    ),
            by = c("target_gene", "target_gene_name"), relationship = "many-to-many"
            ) %>%
  group_by(sRNA, log2FoldChange_mir, pvalue_mir, treatment_mir) %>%
  group_modify(~{
    if(all(is.na(.x$treatment_mrna))){
      .x <- as.data.frame(matrix(NA, nrow = 2, ncol = ncol(.x), dimnames = list(NULL, colnames(.x))))
      .x$treatment_mrna <- c("Venetoclax 3H vs DMSO 3H", "Venetoclax 12H vs DMSO 12H")
    }
    else if(!("Venetoclax 3H vs DMSO 3H" %in% .x$treatment_mrna)){
      .x <- .x[which(!is.na(.x$treatment_mrna)),]
      df <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(.x), dimnames = list(NULL, colnames(.x))))
      df$treatment_mrna <- "Venetoclax 3H vs DMSO 3H"
      .x <- as.data.frame(rbind(df, .x))
    }
    else if(!("Venetoclax 12H vs DMSO 12H" %in% .x$treatment_mrna)){
      .x <- .x[which(!is.na(.x$treatment_mrna)),]
      df <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(.x), dimnames = list(NULL, colnames(.x))))
      df$treatment_mrna <- "Venetoclax 12H vs DMSO 12H"
      .x <- as.data.frame(rbind(.x, df))
    }
    else{
      .x <- .x[which(!is.na(.x$treatment_mrna)),]
    }
    
    return(.x)
  })


mir_up_mapped_mrna_summary <- mir_up_mapped_mrna %>%
  ungroup() %>%
  group_by(sRNA, log2FoldChange_mir, pvalue_mir, treatment_mir, treatment_mrna) %>%
  summarise(N_mrna = sum(!is.na(regulation_mrna)),
            Nup_mrna = sum(regulation_mrna == "up", na.rm = TRUE),
            Nns_mrna = sum(regulation_mrna == "no", na.rm = TRUE),
            Ndown_mrna = sum(regulation_mrna == "down", na.rm = TRUE)) %>%
  filter(!(treatment_mir == "Venetoclax 12H vs DMSO 12H" & treatment_mrna == "Venetoclax 3H vs DMSO 3H"))


### getting probability of having down reg mRNA --> binomial porbability
prob_mir_mrnadown <- function(set, down, ndown, ngene, method = c("exactly", "at least")){
  method <- match.arg(method)
  
  p <- ndown/ngene
  
  if(method == "exactly"){ # probability of having exactly down down-regulated mRNA
    bino <- choose(set, down)                   # binomial coefficient
    res <- bino*((p**down)*((1-p)**(set-down))) # binomial probability
  }
  else if(method == "at least"){# computing 1 minus sum of probability of having exactly 0 to down down-regulated mRNA
    res <- 1 - sum(sapply(0:(down-1), 
                          function(x) choose(set, x)*((p**x)*((1-p)**(set-x)))))
  }
  
  return(res)
}

# test function with miRNA novel_504
prob_mir_mrnadown(mir_up_mapped_mrna_summary$N_mrna[mir_up_mapped_mrna_summary$sRNA == "novel_504" &
                                                      mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                      mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H"],
                  mir_up_mapped_mrna_summary$Ndown_mrna[mir_up_mapped_mrna_summary$sRNA == "novel_504" &
                                                          mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                          mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H"],
                  sum(v12h_dmso12h_diff_mrna$regulation == "down"), nrow(v12h_dmso12h_diff_mrna),
                  "exactly") # 0.01018819

prob_mir_mrnadown(mir_up_mapped_mrna_summary$N_mrna[mir_up_mapped_mrna_summary$sRNA == "novel_504" &
                                                      mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                      mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H"],
                  mir_up_mapped_mrna_summary$Ndown_mrna[mir_up_mapped_mrna_summary$sRNA == "novel_504" &
                                                          mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                          mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H"],
                  sum(v12h_dmso12h_diff_mrna$regulation == "down"), nrow(v12h_dmso12h_diff_mrna),
                  "at least") # 0.02196501
# so quite unlikely to have so many down-regulated mRNA from the target of novel_504
# meaning it must have a real influence on mRNA level


## Compute proba to all
mir_up_mapped_mrna_summary <- mir_up_mapped_mrna_summary %>%
  group_by(sRNA, log2FoldChange_mir, pvalue_mir, treatment_mir, treatment_mrna) %>%
  mutate(p_exact = ifelse(N_mrna == 0, NA,
                          ifelse(treatment_mrna == "Venetoclax 3H vs DMSO 3H",
                                 prob_mir_mrnadown(N_mrna, Ndown_mrna,
                                                   sum(v3h_dmso3h_diff_mrna$regulation == "down"),
                                                   nrow(v3h_dmso3h_diff_mrna), "exactly"),
                                 prob_mir_mrnadown(N_mrna, Ndown_mrna,
                                                   sum(v12h_dmso12h_diff_mrna$regulation == "down"),
                                                   nrow(v12h_dmso12h_diff_mrna), "exactly")
                                 )
                          ),
  p_atleast = ifelse(N_mrna == 0, NA,
                     ifelse(treatment_mrna == "Venetoclax 3H vs DMSO 3H",
                            prob_mir_mrnadown(N_mrna, Ndown_mrna,
                                              sum(v3h_dmso3h_diff_mrna$regulation == "down"),
                                              nrow(v3h_dmso3h_diff_mrna), "at least"),
                            prob_mir_mrnadown(N_mrna, Ndown_mrna,
                                              sum(v12h_dmso12h_diff_mrna$regulation == "down"),
                                              nrow(v12h_dmso12h_diff_mrna), "at least")
                            )
                     )
  )

### plot volcano plot from target mRNA of most influential miRNA
## miRNA 3h and mRNA 12h (considering mRNA decay and longer half life of miRNA, miRNA 3h should have most visible influence on mRNA 12h)
mir3h_influence_mrna12h <- mir_up_mapped_mrna[which(mir_up_mapped_mrna$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                      mir_up_mapped_mrna$treatment_mrna == "Venetoclax 12H vs DMSO 12H" &
                                                      !is.na(match(mir_up_mapped_mrna$sRNA,
                                                                   mir_up_mapped_mrna_summary$sRNA[which(mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                                                                           mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H" &
                                                                                                           mir_up_mapped_mrna_summary$p_exact <= 0.05 &
                                                                                                           mir_up_mapped_mrna_summary$p_atleast <= 0.05)]
                                                                   )
                                                             )
                                                    ),]

mir3h_influence_mrna12h_pannot <- mir_up_mapped_mrna_summary[which(mir_up_mapped_mrna_summary$treatment_mir == "Venetoclax 3H vs DMSO 3H" &
                                                                     mir_up_mapped_mrna_summary$treatment_mrna == "Venetoclax 12H vs DMSO 12H" &
                                                                     mir_up_mapped_mrna_summary$p_exact <= 0.05 &
                                                                     mir_up_mapped_mrna_summary$p_atleast <= 0.05),] %>%
  mutate(x = -0.5, y = 12.5, regulation_mrna = NA,
         lab = paste0("- Mapped mRNA targets: ", N_mrna,
                      "\n- mRNA targets down-regulated: ", Ndown_mrna,
                      "\n- Probability of having at least ", Ndown_mrna, "\ndown-regulated mRNA: ", round(p_atleast, 4)))

# plot mRNA volcano plots from top miRNA
g <- ggplot(mir3h_influence_mrna12h, aes(log2FoldChange_mrna, -log10(pvalue_mrna),
                                         color = factor(regulation_mrna, levels = c("up", "down", "no"))
                                         )) +
  geom_point(size = 2.5, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = mir3h_influence_mrna12h_pannot,
            aes(x, y, label = lab), hjust = 0, fontface = "bold", size = 6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggrepel::geom_label_repel(data = mir3h_influence_mrna12h[-log10(mir3h_influence_mrna12h$pvalue_mrna)*mir3h_influence_mrna12h$log2FoldChange_mrna <= -5 |
                                                             (mir3h_influence_mrna12h$pvalue_mrna <= 0.005 & mir3h_influence_mrna12h$log2FoldChange_mrna <= -0.75) |
                                                             -log10(mir3h_influence_mrna12h$pvalue_mrna) > 5,],
                            aes(log2FoldChange_mrna, -log10(pvalue_mrna), label = target_gene_name,
                                color = factor(regulation_mrna, levels = c("up", "down", "no"))),
                            min.segment.length = 0, show.legend = FALSE, size = 7.5, max.overlaps = 100) +
  facet_wrap(~sRNA) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "no" = "grey30"), na.value = "black") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold-change", y = "-log10(p-value)",
       title = "mRNA Venetoclax 12H vs DMSO 12H\ntargets from the most influential miRNA up-regulated at 3H",
       subtitle = "p-value cutoff: 0.05; log2(fold-change) cutoffs: |1|") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
        plot.subtitle = element_text(face = "bold", size = 22),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(size = 26),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", color = "black", size = 26))

ggsave("Figures/mir3h_influence_mrna12h - mRNA_volcano.png", plot = g,
       device = "png", dpi = 300, height = 12, width = 24)

rm(g)

#### GSEA graph from mRNA targets from three most influential miRNA (Fig 4D)
d <- godata('org.Hs.eg.db', ont="BP")

## performing GSEA 
mir3h_influence_mrna12h_gsea <- mir3h_influence_mrna12h %>% ungroup() %>%
  mutate(score = -log10(pvalue_mrna)*sign(log2FoldChange_mrna)) %>%
  dplyr::select(target_gene, score) %>% unique()
mir3h_influence_mrna12h_gsea <- rlang::set_names(mir3h_influence_mrna12h_gsea$score,
                                          mir3h_influence_mrna12h_gsea$target_gene)
mir3h_influence_mrna12h_gsea <- sort(mir3h_influence_mrna12h_gsea, decreasing = TRUE)
mir3h_influence_mrna12h_gsea <- gseGO(mir3h_influence_mrna12h_gsea,
                                      OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL",
                                      pvalueCutoff = 0.01, minGSSize = 40, maxGSSize = 400, eps = 1e-16,
                                      pAdjustMethod = "BH")


mir3h_influence_mrna12h_gsea_paired <- pairwise_termsim(mir3h_influence_mrna12h_gsea, method = "Lin", semData = d)


mir3h_influence_mrna12h_gsea_hc <- hclust(as.dist(abs(mir3h_influence_mrna12h_gsea_paired@termsim-1)),
                                          method = "average")
# shorter names
mir3h_influence_mrna12h_gsea_hc$labels <- sapply(mir3h_influence_mrna12h_gsea_hc$labels,
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

## prepare graph data
relations <- Reduce(rbind,
                    sapply(unique(mir3h_influence_mrna12h_gsea_paired@result$Description),
                           function(n1){
                             Reduce(rbind,
                                    sapply(unique(mir3h_influence_mrna12h_gsea_paired@result$Description),
                                           function(n2){
                                             if(which(colnames(mir3h_influence_mrna12h_gsea_paired@termsim) == n1) <
                                                which(colnames(mir3h_influence_mrna12h_gsea_paired@termsim) == n2))
                                               data.frame(from = n1, to = n2,
                                                          distance = mir3h_influence_mrna12h_gsea_paired@termsim[n1,n2])
                                           }, simplify = FALSE)
                             )
                           }, simplify = FALSE)
                    ) %>%
  dplyr::filter(distance >= 0.5 & from != to)

relations$from <- sapply(relations$from,
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

relations$to <- sapply(relations$to,
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


names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 1)) # immune response
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 2)) # defense response
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 3)) # translation
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 4)) # cell adhesion
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 5)) # cell activation & differentiation
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 6)) # cytokine production
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 7)) # cytokine response
names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 8)) # "monocarboxylic acid metabolic process"

## adding summary names
relations <- as.data.frame(rbind(relations,
                                 rbind(data.frame(from = "Immune response",
                                                  to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 1)),
                                                  distance = 1),
                                       rbind(data.frame(from = "Defense response",
                                                        to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 2)),
                                                        distance = 1),
                                             rbind(data.frame(from = "Translation",
                                                              to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 3)),
                                                              distance = 1),
                                                   rbind(data.frame(from = "Cell adhesion",
                                                                    to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 4)),
                                                                    distance = 1),
                                                         rbind(data.frame(from = "Cell activation & differentiation",
                                                                          to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 5)),
                                                                          distance = 1),
                                                               rbind(data.frame(from = "Cytokine production",
                                                                                to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 6)),
                                                                                distance = 1),
                                                                     data.frame(from = "Cytokine response",
                                                                                to = names(which(cutree(mir3h_influence_mrna12h_gsea_hc, 8) == 7)),
                                                                                distance = 1)
                                                               )
                                                         )
                                                   )
                                             )
                                       )
                                 )
                                 )
                           )

relations <- as.data.frame(rbind(relations,
                                 data.frame(from = c("Immune response", "Cell adhesion", "Cytokine response",
                                                     "Immune response", "Immune response", "Defense response"),
                                            to = c("Defense response", "Cell activation & differentiation",
                                                   "Cytokine production", "Cytokine response", "Cytokine production",
                                                   "Cytokine production"),
                                            distance = c(1,1,1,0.9,0.9,0.9)
                                            )
                                 )
                           )


vertex <- mir3h_influence_mrna12h_gsea_paired@result %>%
  dplyr::select(Description, enrichmentScore, NES, p.adjust, core_enrichment) %>%
  group_by(Description, enrichmentScore, NES, p.adjust) %>%
  summarise(GeneCount = length(strsplit(core_enrichment, "/")[[1]]))

vertex$Description <- sapply(vertex$Description,
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
vertex$labsize <- 1

vertex <- as.data.frame(rbind(vertex,
                              data.frame(Description = c("Immune response", "Defense response",
                                                         "Cell activation & differentiation", "Cell adhesion",
                                                         "Cytokine response", "Cytokine production", "Translation"),
                                         GeneCount = sapply(c("Immune response", "Defense response",
                                                              "Cell activation & differentiation", "Cell adhesion",
                                                              "Cytokine response", "Cytokine production", "Translation"),
                                                            function(x) length(unique(unlist(strsplit(mir3h_influence_mrna12h_gsea_paired@result$core_enrichment[which(!is.na(match(mir3h_influence_mrna12h_gsea_paired@result$Description,
                                                                                                                                                                                    gsub("\n", " ", relations$to[relations$from == x])
                                                                                                                                                                                    )))
                                                                                                                                                                 ], "/")))
                                                                               ),
                                                            USE.NAMES = FALSE
                                                            ),
                                         enrichmentScore = sapply(c("Immune response", "Defense response",
                                                                    "Cell activation & differentiation", "Cell adhesion",
                                                                    "Cytokine response", "Cytokine production", "Translation"),
                                                                  function(x) median(vertex$enrichmentScore[which(!is.na(match(vertex$Description,
                                                                                                                               gsub("\n", " ", relations$to[relations$from == x])
                                                                                                                               )))
                                                                                                            ]),
                                                                  USE.NAMES = FALSE
                                                                  ),
                                         NES = sapply(c("Immune response", "Defense response",
                                                        "Cell activation & differentiation", "Cell adhesion",
                                                        "Cytokine response", "Cytokine production", "Translation"),
                                                      function(x) median(vertex$NES[which(!is.na(match(vertex$Description,
                                                                                                       gsub("\n", " ", relations$to[relations$from == x])
                                                                                                       )))
                                                                                    ]),
                                                      USE.NAMES = FALSE
                                                      ),
                                         p.adjust = sapply(c("Immune response", "Defense response",
                                                             "Cell activation & differentiation", "Cell adhesion",
                                                             "Cytokine response", "Cytokine production", "Translation"),
                                                           function(x) metap::logitp(vertex$p.adjust[which(!is.na(match(vertex$Description,
                                                                                                                        gsub("\n", " ", relations$to[relations$from == x])
                                                                                                                        )))
                                                                                                     ])$p,
                                                           USE.NAMES = FALSE
                                                           ),
                                         labsize = 100
                                         )
                              )
                        )


## plotting graph
g <- igraph::graph_from_data_frame(relations, directed = FALSE, vertices = vertex)

p <- ggraph(g, layout = "stress", weights = distance) +
  geom_edge_link(show.legend = FALSE,
                 aes(alpha = distance, color = distance, width = distance)) +
  geom_node_point(aes(color = p.adjust, fill = NES, size = GeneCount),
                  shape = 21, stroke = 1.8) +
  geom_node_label(aes(x = x, y = y, label = name, size = labsize),
                  color = "black", show.legend = FALSE, fontface = "bold",
                  min.segment.length = 0, max.overlaps = 100, repel = TRUE) +
  scale_size_identity() +
  scale_fill_gradient(high = "white", low = "blue", name = "NES", limits = c(0,-2)) +
  scale_color_gradient(name = "p.adjust", transform = "log", high = "#E6A540", low = "#00890B") +
  scale_size_continuous(name = "Gene count", limits = c(0,400),
                        breaks = c(5,10,20,60,seq(100,200, 50)),
                        range = c(2,16)) +
  scale_edge_alpha_continuous(range = c(0.2,0.8)) +
  scale_edge_width_continuous(range = c(0.5,2)) +
  scale_edge_color_gradient(low = "white", high = "red") +
  labs(title = "GSEA graph terms from novel_504, miR-1039-3p and miR-1039-5p targets",
       subtitle = "adj.p-value cutoff: 0.01 - minimum genes: 10") +
  theme_graph(base_family="sans") +
  theme(legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14),
        plot.title = element_text(face = "bold", color = "black", size = 20),
        plot.subtitle = element_text(size = 18))

p@layers[[3]]$aes_params$colour <- c(rep("grey40", length(mir3h_influence_mrna12h_gsea_hc$labels)), 
                                     rep("#a7c957", 7))
p@layers[[3]]$aes_params$fill <- c(rep("white", length(mir3h_influence_mrna12h_gsea_hc$labels)), 
                                   rep("#f2e8cf", 7))

ggsave("Figures/mir3h_influence_mrna12h - GSEA_graph.png", plot = p, 
       dpi = 300, width = 18, heigh = 10)

# data cleaning
rm(d, mir3h_influence_mrna12h_gsea, mir3h_influence_mrna12h_gsea_paired, 
   mir3h_influence_mrna12h_gsea_hc, relations, vertex, g, p)


#### Volcanos from most influential miRNAs (SFig 8D)
## 3h
g <- ggplot(v3h_dmso3h_diff_mir, aes(log2FoldChange, -log10(pval),
                                 color = factor(regulation, levels = c("up", "down", "no"))
                                 )) +
  geom_point(size = 1.5) +
  geom_point(data = v3h_dmso3h_diff_mir %>%
               dplyr::filter(sRNA %in% c("novel_504", "hsa-miR-10394-3p", "hsa-miR-10394-5p")),
             aes(log2FoldChange, -log10(pval),
                 color = factor(regulation, levels = c("up", "down", "no"))),
             shape = 24, size = 3, stroke = 1.4, fill = "orange",
             show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggrepel::geom_label_repel(data = v3h_dmso3h_diff_mir %>%
                              dplyr::filter(sRNA %in% c("novel_504", "hsa-miR-10394-3p", "hsa-miR-10394-5p")),
                            aes(log2FoldChange, -log10(pval), label = sRNA,
                                color = factor(regulation, levels = c("up", "down", "no"))),
                            min.segment.length = 0, show.legend = FALSE, size = 7.5, max.overlaps = 50) +
  ylim(c(0,9)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-6,6),
                     labels = paste(round(c(-6, -4, -2, 0, 2, 4, 6), 1))) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "no" = "grey30"),
                     labels = c("up" = paste0("up-regulated (n=", sum(v3h_dmso3h_diff_mir$regulation == "up"), ")"),
                                "down" = paste0("down-regulated (n=", sum(v3h_dmso3h_diff_mir$regulation == "down"), ")"),
                                "no" = paste0("non-significant (n=", sum(v3h_dmso3h_diff_mir$regulation == "no"), ")")),
                     name = "p-value <= 0.05\n|log2(fold-change)| >= 1") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold-change", y = "-log10(p-value)",
       title = "Venetoclax 3H vs DMSO 3H") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(size = 26),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(size = 26))

ggsave("Figures/volcano_V3HvsDMSO3H_miR_MostInfluential.png", plot = g,
       device = "png", dpi = 300, height = 12, width = 16)

## Venetoclax 12h vs Venetoclax 3h
v12h_v3h_diff_mir <- readr::read_tsv("Data/RNAseq/sRNA/V12HvsV3H.Differential_analysis_results.xls")
v12h_v3h_diff_mir$regulation <- ifelse(v12h_v3h_diff_mir$pval <= 0.05,
                                           ifelse(v12h_v3h_diff_mir$log2FoldChange > 1, "up",
                                                  ifelse(v12h_v3h_diff_mir$log2FoldChange < -1, "down", "no")),
                                           "no")
v12h_v3h_diff_mir$regulation[which(is.na(v12h_v3h_diff_mir$regulation))] <- "no"

g <- ggplot(v12h_v3h_diff_mir, aes(log2FoldChange, -log10(pval),
                                   color = factor(regulation, levels = c("up", "down", "no"))
                                   )) +
  geom_point(size = 1.5) +
  geom_point(data = v12h_v3h_diff_mir %>%
               dplyr::filter(sRNA %in% c("novel_504", "hsa-miR-10394-3p", "hsa-miR-10394-5p")),
             aes(log2FoldChange, -log10(pval),
                 color = factor(regulation, levels = c("up", "down", "no"))),
             shape = 24, size = 3, stroke = 1.4, fill = "orange",
             show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggrepel::geom_label_repel(data = v12h_v3h_diff_mir %>%
                              dplyr::filter(sRNA %in% c("novel_504", "hsa-miR-10394-3p", "hsa-miR-10394-5p")),
                            aes(log2FoldChange, -log10(pval), label = sRNA,
                                color = factor(regulation, levels = c("up", "down", "no"))),
                            min.segment.length = 0, show.legend = FALSE, size = 7.5, max.overlaps = 50) +
  ylim(c(0,9)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-6,6),
                     labels = paste(round(c(-6, -4, -2, 0, 2, 4, 6), 1))) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "no" = "grey30"),
                     labels = c("up" = paste0("up-regulated (n=", sum(v12h_v3h_diff_mir$regulation == "up"), ")"),
                                "down" = paste0("down-regulated (n=", sum(v12h_v3h_diff_mir$regulation == "down"), ")"),
                                "no" = paste0("non-significant (n=", sum(v12h_v3h_diff_mir$regulation == "no"), ")")),
                     name = "p-value <= 0.05\n|log2(fold-change)| >= 1") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold-change", y = "-log10(p-value)",
       title = "Venetoclax 12H vs Venetoclax 3H") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(size = 26),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(size = 26))

ggsave("Figures/volcano_V12HvsV3H_miR_MostInfluential.png", plot = g,
       device = "png", dpi = 300, height = 12, width = 16)

rm(g)

#### Venn SFig 8B - C
### mRNA target Venn diagrams
venn <- VennDiagram::venn.diagram(sapply(unique(mir3h_influence_mrna12h$sRNA),
                                         function(x) mir3h_influence_mrna12h$target_gene[which(mir3h_influence_mrna12h$sRNA == x)],
                                         simplify = FALSE),
                                  filename = NULL, disable.logging = TRUE, margin = 0.5,
                                  print.mode = "percent",
                                  fill = c("red", "blue", "green"), lwd = 2, cex = 3,
                                  main = "Venn Diagram of predicted mRNA targets from\nmiRNA most likely to influence mRNA level at 12h",
                                  main.fontface = "bold", main.cex = 3, main.pos = c(0.5,0.8),
                                  #cat.pos = c(220, 180),
                                  cat.dist = c(0.08,0.08, 0.04),
                                  fontface = "bold", fontfamiliy = "sans",
                                  cat.cex = 2.4, cat.fontface = "bold",
                                  cat.fontfamily = "sans")
png("Figures/venn_3miRNA_target_translation_all.png",
    width=12, height=12, units = "in", res=300)
grid::grid.newpage()
grid::grid.draw(venn)
dev.off()


venn <- VennDiagram::venn.diagram(sapply(unique(mir3h_influence_mrna12h$sRNA),
                                         function(x) mir3h_influence_mrna12h$target_gene[which(mir3h_influence_mrna12h$sRNA == x &
                                                                                                 mir3h_influence_mrna12h$regulation_mrna == "down")],
                                         simplify = FALSE),
                                  filename = NULL, disable.logging = TRUE, margin = 0.5,
                                  print.mode = "raw",
                                  fill = c("red", "blue", "green"), lwd = 2, cex = 3,
                                  main = "Venn Diagram of predicted mRNA targets\ndown-regulated from\nmiRNA most likely to influence mRNA level at 12h",
                                  main.fontface = "bold", main.cex = 3, main.pos = c(0.5,0.8),
                                  #cat.pos = c(220, 180),
                                  cat.dist = c(0.08,0.08, 0.04),
                                  fontface = "bold", fontfamiliy = "sans",
                                  cat.cex = 2.4, cat.fontface = "bold",
                                  cat.fontfamily = "sans")
png("Figures/venn_3miRNA_target_translation_DownOnly.png",
    width=12, height=12, units = "in", res=300)
grid::grid.newpage()
grid::grid.draw(venn)
dev.off()



#### saving data for supplementary material
openxlsx::write.xlsx(list("V3HvsDMSO3H" = v3h_dmso3h_diff_mrna,
                          "V12HvsDMSO12H" = v12h_dmso12h_diff_mrna),
                     "Results/DifferentialAnalysis_mRNA_VenetoclaxKasumi-1.xlsx")

openxlsx::write.xlsx(list("V3HvsDMSO3H" = v3h_dmso3h_diff_mir,
                          "V12HvsDMSO12H" = v12h_dmso12h_diff_mir,
                          "V12HvsV3H" = v12h_v3h_diff_mir),
                     "Results/DifferentialAnalysis_smallRNA_VenetoclaxKasumi-1.xlsx")
openxlsx::write.xlsx(full_join(mir_up_mapped_mrna_summary,
                               mir_up_mapped_mrna,
                               by = intersect(colnames(mir_up_mapped_mrna),
                                              colnames(mir_up_mapped_mrna_summary))
                               ) %>%
                       dplyr::rename(binomial_p_exact = p_exact, binomial_p_atleast = p_atleast),
                     "Results/miRNA_up-reg_PredTarget_BinoProb.xlsx" 
                     )


