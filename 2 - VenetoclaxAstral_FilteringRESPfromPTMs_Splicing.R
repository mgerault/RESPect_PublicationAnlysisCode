library(IMPRINTS.CETSA.app)

vene_astral_catasresp <- as.data.frame(rbind(vene_astral_resp3h[which(vene_astral_resp3h$category == "RESP"),],
                                             vene_astral_resp12h[which(vene_astral_resp12h$category == "RESP"),]
                                             )
                                       )
length(unique(vene_astral_catasresp$id)) # 100 proteins categorized as RESP

# removing Nvalue and Npep columns
vene_astral_catasresp <- vene_astral_catasresp[,-grep("^N", colnames(vene_astral_catasresp))]

## which could be isoform ?
vene_astral_catasresp$Gene[vene_astral_catasresp$id %in% 
                             unique(c(vene_astral_isoform3h$accession,
                                      vene_astral_isoform12h$accession)
                                    )] # 5
# ATRX, CIAPIN1 and HNRNPA3 at 3h (HNRNPA3 was categorized as Unconfident at 3h)
# HNRNPA3, PPP1R8 and ZMYND8 at 12h (HNRNPA3 was categorized as RESP at 12h)

# from plots, ATRX and HNRNPA3 cannot be due to alternative splicing at 3h
# however, CIAPIN1 is more likely to be due to its second isoform being more expressed

# from plots, PPP1R8, HNRNPA3 and ZMYND8 cannot be due to alternative splicing at 12h

# will remove CIAPIN1 from the cleavage list

## which has shifting peptides due to PTMs ?
length(unique(vene_astral_catasresp$Gene[vene_astral_catasresp$id %in% 
                                           unique(c(vene_astral_ptm3h$id,
                                                    vene_astral_ptm12h$id))
                                         ]
              )) # 90
 
### Do same clustering as done for isoform (up, down, ambiguous)
# check proportion of mapped PTMs in down
# if proportion high, check up cluster for proportion of quantified modif and/or value compare to down
# if down has greater FC than up, PTMs are more likely than cleavage
vene_astral_catasresp_PTMdata <- imprints_join_peptides(list(vene_astral_pepdiff3h[which(!is.na(match(vene_astral_pepdiff3h$Master.Protein.Accessions, 
                                                                                                   unique(vene_astral_catasresp$id)))),
                                                                                -grep("Vehicle", colnames(vene_astral_pepdiff3h))], 
                                                          vene_astral_pepdiff12h[which(!is.na(match(vene_astral_pepdiff12h$Master.Protein.Accessions, 
                                                                                                   unique(vene_astral_catasresp$id)))),
                                                                                -grep("Vehicle", colnames(vene_astral_pepdiff12h))]),
                                                     "exact")
vene_astral_catasresp_PTMdata$Annotated.Sequence <- NULL
vene_astral_catasresp_PTMdata$Gene <- sapply(vene_astral_catasresp_PTMdata$description, 
                                          IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE)
vene_astral_catasresp_PTMdata$description <- sapply(vene_astral_catasresp_PTMdata$description, 
                                                 IMPRINTS.CETSA.app:::getProteinName, USE.NAMES = FALSE)
colnames(vene_astral_catasresp_PTMdata)[1] <- "id"

vene_astral_catasresp_PTMdata <- vene_astral_catasresp_PTMdata %>%
  tidyr::gather("treatment", "value", -id, -Positions.in.Master.Proteins,
                -description, -Gene, -Modifications) %>%
  tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(id, description, Gene, Positions.in.Master.Proteins, Modifications, temp, treatment) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% dplyr::group_by(id, description, Gene, Positions.in.Master.Proteins, Modifications, treatment) %>%
  # compute greatest FC from all temperatures
  dplyr::reframe(maxFC = mean_value[which.max(abs(mean_value))]) %>%
  # cluster peptides as in imprints_categorize_peptides but not on absolute value
  dplyr::ungroup() %>% dplyr::group_by(id, description, Gene, treatment) %>%
  dplyr::group_modify(~ {
    h <- data.frame(comp = cutree(hclust(dist(.x$maxFC), method = "complete"), 2),
                    ward = cutree(hclust(dist(.x$maxFC), method = "ward.D"), 2),
                    ward2 = cutree(hclust(dist(.x$maxFC), method = "ward.D2"), 2),
                    ave = cutree(hclust(dist(.x$maxFC), method = "average"), 2))
    h_name <- apply(h, 2,
                    function(y){
                      sapply(1:2, function(z) mean(.x$maxFC[y == z]))
                    })
    if(!all(apply(h_name, 2, which.max) == 1)){
      h[,which(apply(h_name, 2, which.max) == 2)] <-
        apply(h[,which(apply(h_name, 2, which.max) == 2),drop=FALSE], 2,
              function(y){
                y[y == 1] <- 3
                y[y == 2] <- 1
                y[y == 3] <- 2;
                y
              })
    }
    h <- apply(h, 1, mean)
    
    if(any(h == 1.5)){ # ambiguous peptides will be in cluster named 3
      p_ambi <- which(h == 1.5)
      h[p_ambi] <- 3
    }
    h <- round(h)
    
    .x$cluster <- h
    
    return(.x)
  })

vene_astral_catasresp_PTMdata$Modifications <- mapply(function(m, p){
  m <- strsplit(m, "\\]; ")[[1]]
  m <- m[-grep("TMTpro", m)]
  if(length(m)){
    p <- sub(";.*", "", p)
    p <- gsub(".* \\[|-.*", "", p)
    p <- as.numeric(p) - 1
    
    m_change <- gsub(".*\\[|\\]", "", m)
    m_change <- unlist(strsplit(m_change, "; "))
    m_change <- sapply(m_change,
                       function(y){
                         yp <- sub("^.", "", y)
                         y <- sub(yp, "", y)
                         if(nchar(yp)){
                           yp <- as.numeric(yp) + p
                         }
                         else{
                           yp <- p
                         };
                         paste0(y, yp)
                       })
    
    m <- paste(m, collapse = "]; ")
    if(!grepl("\\]$", m)){
      m <- paste0(m, "]")
    }
    
    for(i in names(m_change)){
      m <- sub(paste0("(?<=(\\[| ))", i), m_change[[i]], m, perl = TRUE)
    }
  }
  else{
    m <- NA
  };
  m
}, vene_astral_catasresp_PTMdata$Modifications, 
vene_astral_catasresp_PTMdata$Positions.in.Master.Proteins, 
USE.NAMES = FALSE)

vene_astral_catasresp_PTMdata$Modifications.Position <- sapply(vene_astral_catasresp_PTMdata$Modifications,
                                      function(x){
                                        if(!is.na(x)){
                                          x <- strsplit(x, "\\]; ")[[1]]
                                          x <- gsub(".* \\[|\\]", "", x)
                                          x <- paste(x, collapse = "; ")
                                          x <- strsplit(x, "; ")[[1]]
                                          x <- x[order(as.numeric(sub("^.", "", x)))]
                                          x <- paste(x, collapse = "; ")
                                        };
                                        x
                                      }, USE.NAMES = FALSE)
vene_astral_catasresp_PTMdata$Modifications <- sapply(vene_astral_catasresp_PTMdata$Modifications,
                             function(x){
                               if(!is.na(x)){
                                 x <- strsplit(x, "\\]; ")[[1]]
                                 x <- sub(" .*", "", x)
                                 x <- sub("^\\d{1}x", "", x)
                                 x <- x[order(x)]
                                 x <- paste(x, collapse = "; ")
                               };
                               x
                             }, USE.NAMES = FALSE)

vene_astral_catasresp_PTMdata <- left_join(vene_astral_catasresp_PTMdata, 
                                        as.data.frame(rbind(vene_astral_ptm3h[,c(1:7,14:17)],
                                                            vene_astral_ptm12h[,c(1:7,14:17)])),
                              by = c("id", "Gene", "description", "Positions.in.Master.Proteins",
                                     "Modifications", "Modifications.Position", "treatment"),
                              relationship = "many-to-many")
vene_astral_catasresp_PTMdata <- left_join(vene_astral_catasresp, vene_astral_catasresp_PTMdata,
                              by = c("id", "Gene", "description", "treatment"),
                              relationship = "one-to-many")

modif_abbrev <- list("ACETYLATION" = "ac", "SUCCINYLATION" = "sc", "PHOSPHORYLATION" = "p",
                     "UBIQUITYLATION" = "ub", "SUMOYLATION" = "sm",
                     "MONO-METHYLATION" = "m1", "DI-METHYLATION" = "m2",
                     "CASPASE CLEAVAGE" = "ca")
vene_astral_catasresp_PTMdata$PTM <- mapply(function(x, a, p){
                                 if(is.na(x)){
                                   NA
                                 }
                                 else{
                                   paste0(a, p, "-", modif_abbrev[[x]])
                                 }
                               }, vene_astral_catasresp_PTMdata$PTM, vene_astral_catasresp_PTMdata$PTM_aa, vene_astral_catasresp_PTMdata$PTM_position, 
                               USE.NAMES = FALSE)

vene_astral_catasresp_PTMdata$PTM_aa <- NULL
vene_astral_catasresp_PTMdata$PTM_position <- NULL
vene_astral_catasresp_PTMdata$Modifications <- NULL
colnames(vene_astral_catasresp_PTMdata)[13] <- "Modifications"

vene_astral_catasresp_PTMdata <- vene_astral_catasresp_PTMdata %>% dplyr::ungroup() %>%
  dplyr::group_by(id, Gene, description, treatment, combined_pvalue, RESP_score,
                  cleaved_site, category, details,
                  Positions.in.Master.Proteins, maxFC, cluster) %>%
  summarise(Modifications = paste(na.omit(unique(Modifications)), collapse = "; "),
            PTM = paste(na.omit(PTM), collapse = "; "),
            nbref = paste(na.omit(nbref), collapse = "; "))

vene_astral_catasresp_PTMdata$Modifications[which(vene_astral_catasresp_PTMdata$Modifications == "")] <- NA
vene_astral_catasresp_PTMdata$PTM[which(vene_astral_catasresp_PTMdata$PTM == "")] <- NA
vene_astral_catasresp_PTMdata$nbref[which(vene_astral_catasresp_PTMdata$nbref == "")] <- NA

vene_astral_catasresp_PTMprop <- vene_astral_catasresp_PTMdata %>% dplyr::ungroup() %>%
  dplyr::group_by(id, Gene, description, treatment, combined_pvalue, RESP_score,
                  cleaved_site, category, details, cluster) %>%
  group_modify(~ {
    m <- mean(.x$maxFC)
    ds <- 1/(1 + (.x$maxFC - m)**2)
    prop <- ds/sum(ds)
    
    ptm <- which(is.na(.x$Modifications) & !is.na(.x$PTM))
    modif <- which(!is.na(.x$Modifications) & is.na(.x$PTM))
    ptm_modif <- which(!is.na(.x$Modifications) & !is.na(.x$PTM))
    none <- which(is.na(.x$Modifications) & is.na(.x$PTM))
    
    ptm <- ifelse(length(ptm), sum(prop[ptm]), 0)
    modif <- ifelse(length(modif), sum(prop[modif]), 0)
    ptm_modif <- ifelse(length(ptm_modif), sum(prop[ptm_modif]), 0)
    none <- ifelse(length(none), sum(prop[none]), 0)
    
    df <- data.frame(npep = nrow(.x), PTM = ptm, Modified = modif, PTM_Modified = ptm_modif, None = none)
    
    return(df)
  })

# From the data, definitely not cleavage but PTM are:
vene_astral_catasresp_PTMprop %>% dplyr::ungroup() %>%
  dplyr::group_by(id, Gene, description, treatment, combined_pvalue, RESP_score,
                  cleaved_site,  category, details) %>%
  filter(PTM[cluster == 2] + PTM_Modified[cluster == 2] >= 0.5 & 
           None[cluster == 1] < 0.75 & 
           length(unique(cluster)) == 2) %>% 
  dplyr::ungroup() %>% dplyr::select(Gene) %>% unique() %>% unlist() %>% unname() %>% sort()
# ACADS, CASP6, TPP2

# From the data, might be cleavage but Modification are:
vene_astral_catasresp_PTMprop %>% dplyr::ungroup() %>%
  dplyr::group_by(id, Gene, description, treatment, combined_pvalue, RESP_score,
                  cleaved_site,  category, details) %>%
  filter(Modified[cluster == 1] + PTM_Modified[cluster == 1] >= 0.5 & 
           Modified[cluster == 2] + PTM_Modified[cluster == 2] <= 0.5 & 
           length(unique(cluster)) == 2) %>% 
  dplyr::ungroup() %>% dplyr::select(Gene) %>% unique() %>% unlist() %>% unname() %>% sort()
# FGD3

# Most likely cleaved protein without any reported PTMs:
sort(unique(vene_astral_catasresp$Gene)[!(unique(vene_astral_catasresp$Gene) %in% 
                                                                      unique(c(vene_astral_ptm3h$Gene,
                                                                               vene_astral_ptm12h$Gene))
                                                                    )]) # 10
# "CDKN2AIP" "LETM1" "NAT10"   "NOP14" "NOP2"  "PPP1R8"  "PPP1R9B"  "PRDX4" "SUGP1"  "USP36"


### Filter out
vene_astral_catasresp_false <- c("CIAPIN1", # due to alternative splicing
                                 "ACADS", "CASP6", "TPP2", # due to PTMs (oxidation, phospho, ub, ac)
                                 "FGD3" # due to PTMs (heavily carbamidomethylated)
                                 )

vene_astral_filteredresp <- vene_astral_catasresp[which(is.na(match(vene_astral_catasresp$Gene, 
                                                                    vene_astral_catasresp_false))),]
length(unique(vene_astral_filteredresp$id)) 
# 95 RESP which cannot be explained solely by PTMs and/or alternative splicing
