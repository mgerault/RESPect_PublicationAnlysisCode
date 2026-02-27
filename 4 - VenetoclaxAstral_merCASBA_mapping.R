library(IMPRINTS.CETSA.app)

### load merCASBA
mercasba <- openxlsx::read.xlsx("Data/mapped_merCASBA.xlsx")

### Function to get sequence from computed cleavage site
sequence_cleavedsite <- function(data, fasta = "", margin = 10){
  if(!grepl("\\.fasta$", fasta) | !file.exists(fasta)){
    message("Error: fasta must be a valid path to a FASTA file ! No sequences were added to the data.")
    return(data)
  }
  if(!all(c("id", "Gene", "treatment", "cleaved_site") %in% colnames(data))){
    missing_column <- c("id", "Gene", "treatment", "cleaved_site")[!(c("id", "Gene", "treatment", "cleaved_site") %in% colnames(data))]
    message(paste("Error:", paste(missing_column, collapse = ", "), " column(s) is/are mssing in your data ! Your original data has been returned."))
    return(data)
  }
  
  fasta <- seqinr::read.fasta(fasta, seqtype = "AA")
  fasta <- lapply(fasta, as.character)
  names(fasta) <- sapply(strsplit(names(fasta), "\\|"), "[[", 2)
  
  cleavedsite_to_sequence <- function(protein, cleaved_site, margin = 10){
    protein <- strsplit(protein, "; ")[[1]]
    sequence <- fasta[protein]
    
    for(i in 1:length(protein)){
      sequencei <- sequence[[i]]
      if(length(sequencei)){
        cleaved_sitei <- strsplit(cleaved_site, "; ")[[1]][i]
        cleaved_sitei <- strsplit(cleaved_sitei, "-|~")[[1]]
        cleaved_sitei <- as.numeric(cleaved_sitei)
        cleaved_sitei <- cleaved_sitei[order(cleaved_sitei)]
        cleaved_sitei[1] <- cleaved_sitei[1] - margin
        cleaved_sitei[2] <- cleaved_sitei[2] + margin
        if(cleaved_sitei[1] < 1){
          cleaved_sitei[1] <- 1
        }
        if(cleaved_sitei[2] > length(sequencei)){
          cleaved_sitei[2] <- length(sequencei)
        }
        
        sequence[[i]] <- paste(sequencei[cleaved_sitei[1]:cleaved_sitei[2]], collapse = "")
      }
      else{
        message(paste("Warning:", protein[i], "wasn't found in your FASTA file !"))
        sequence[[i]] <- NA
      }
    }
    sequence <- paste(sequence, collapse = "; ")
    return(sequence)
  }
  
  data <- data %>%
    group_by(id, Gene, treatment, cleaved_site) %>%
    mutate(cleaved_sequence = cleavedsite_to_sequence(id, cleaved_site, margin = margin))
  
  return(data)
}

### mapping filtered RESP results
vene_astral_resp_final <- sequence_cleavedsite(vene_astral_resp_final, 
                                               "Data/20240619_UP000005640_9606.fasta")

## handling protein group
unique(grep("; ", vene_astral_resp_final$id, value = TRUE)) 
# "P38159; Q96E39" --> RBMX
sapply(unique(grep("; ", vene_astral_resp_final$id, value = TRUE)), 
       function(x){
         x <- strsplit(x, "; ")[[1]]
         x %in% mercasba$id
       }) # none are in merCASBA
# so can join datasets safely -->

vene_astral_resp_final <- left_join(vene_astral_resp_final, 
                                    mercasba[,c("id", "Gene", "cleaved.pos", 
                                                "cleaved.sequence", "nb.ref")] %>%
                                      dplyr::rename(cleaved_sequence.merCASBA = cleaved.sequence,
                                                    cleaved_site.merCASBA = cleaved.pos,
                                                    nb_ref.merCASBA = nb.ref), 
                                    by = c("id", "Gene"), 
                                    relationship = "many-to-many")


# adding information
vene_astral_resp_final$in_merCASBA <- !is.na(vene_astral_resp_final$cleaved_site.merCASBA)
vene_astral_resp_final$cleaved_site_in_merCASBA <- mapply(function(x, y){
  if(is.na(y)){
    res <- FALSE
  }
  else{
    x <- strsplit(x, "; ")[[1]]
    res <- sapply(x, function(z){
      z <- strsplit(z, "-|~")[[1]]
      z <- as.numeric(z)
      z <- z[order(z)];
      z[1] <= y & y <= z[2]
    }, USE.NAMES = FALSE)
  };
  any(res)
}, vene_astral_resp_final$cleaved_site, vene_astral_resp_final$cleaved_site.merCASBA,
USE.NAMES = FALSE)



vene_astral_resp_final <- vene_astral_resp_final[,c("id", "Gene", "description", "time_category", "treatment", "combined_pvalue", "RESP_score",
                                                    "significant_shift", "category", "details", "cleaved_site", "cleaved_sequence",
                                                    "in_merCASBA", "cleaved_site_in_merCASBA", 
                                                    "cleaved_site.merCASBA", "cleaved_sequence.merCASBA", "nb_ref.merCASBA")]


length(unique(vene_astral_resp_final$id)) # 95 unique proteins
vene_astral_resp_final %>% group_by(id, Gene) %>%
  summarise(in_merCASBA = any(in_merCASBA),
            site_in_merCASBA = any(cleaved_site_in_merCASBA)) %>%
  ungroup() %>%
  summarise(in_merCASBA = sum(in_merCASBA), # 38 / 95
            site_in_merCASBA = sum(site_in_merCASBA) # 10 / 38
            )


### Adding references from all proteins reported to be cleaved by merCASBA
vene_astral_resp_final <- readxl::read_xls("Data/merCASBA.xls") %>%
  dplyr::rename(Uniprot.Gene = `UniProt ID`) %>%
  dplyr::select(-Pos, -Sequence, -`P6-P4'`) %>%
  full_join(openxlsx::read.xlsx("Data/mapped_merCASBA.xlsx") , 
            by = "Uniprot.Gene", relationship = "many-to-many") %>%
  dplyr::select(-Protein.Name, -cleaved.site, -Uniprot.Gene) %>%
  dplyr::rename(cleaved_site.merCASBA = cleaved.pos,
                cleaved_sequence.merCASBA = cleaved.sequence) %>%
  right_join(vene_astral_resp_final,
             by = c("id", "Gene", "cleaved_site.merCASBA", "cleaved_sequence.merCASBA"),
             relationship = "many-to-many") %>% 
  dplyr::select(-nb.ref)

# checking the only mapped references used by merCASBA
vene_astral_resp_final <- vene_astral_resp_final %>%
  group_by(id, description, Gene, time_category, treatment, combined_pvalue, RESP_score,
           significant_shift, cleaved_site, cleaved_sequence, category, details, in_merCASBA,
           cleaved_site_in_merCASBA, cleaved_site.merCASBA, cleaved_sequence.merCASBA, nb_ref.merCASBA
           ) %>%
  group_modify(~{
    ref <- unique(na.omit(unlist(.x)))
    if(length(ref)){
      ref <- sub("\\(.*", "", ref)
      pmid <- sapply(ref, 
                     function(x){
                       mercasbapmid <- openxlsx::read.xlsx("Data/merCASBA_references.xlsx")
                       mercasbapmid$pmid[grep(x, mercasbapmid$references)]
                     }, USE.NAMES = FALSE)
      
      df <- data.frame(merCASBA_reference = paste(sub(" $", "", ref), collapse = "; "), 
                       merCASBA_PMID = ifelse(any(!is.na(pmid)),
                                              paste(pmid, collapse = "; "),
                                              NA)
                       )
    }
    else
      df <- data.frame(merCASBA_reference = NA, merCASBA_PMID = NA)
    
    return(df)
  })
length(unique(unlist(strsplit(na.omit(vene_astral_resp_final$merCASBA_reference), "; ")))) 
# 18 unique references
length(unique(unlist(strsplit(na.omit(vene_astral_resp_final$merCASBA_PMID), "; ")))) 
# 17 PMID (one reference uncheckable: M. Hoek, personal communication)

vene_astral_resp_final <- vene_astral_resp_final[order(vene_astral_resp_final$Gene),]

# openxlsx::write.xlsx(vene_astral_resp_final, "Results/VenetoclaxAstral_RESP_final_summary.xlsx")
