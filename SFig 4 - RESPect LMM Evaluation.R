library(IMPRINTS.CETSA.app)

#### LMM evaluation --> rewrite imprints_cleaved_peptides function
## see doc with ?imprints_cleaved_peptides

#### Internal functions from IMPRINTS.CETSA.app
### function to compute discrete inflexion point
inflex <- function(mp){
  ### compute cumulative sum
  if(any(is.na(mp)))
    mp[which(is.na(mp))] <- 0
  
  x <- cumsum(mp)
  
  ### compute 'exact' theoretical inflexion point using second derivative
  x_df1 <- diff(x) # get first 'first derivative'
  # compute 'tendency' from 'first derivative', i.e. if increasing/decreasing
  df <- data.frame(x = 1:length(x_df1),
                   y = x_df1)
  res <- lm(y ~ x, data = df)
  res <- res$coeff[2]
  res <- sign(res)
  
  x_df2 <- diff(x_df1) # get second 'first derivative'
  if(res == 1){
    inflex_point_exac <- which.max(x_df2) + 2 # double diff, so add 2 to index point
  }
  else{
    inflex_point_exac <- which.min(x_df2) + 2 # double diff, so add 2 to index point
  }
  
  if(inflex_point_exac - 2 == length(x_df2)){  # if inflex_point is last one, can't conclude that is cleaved
    inflex_point_exac <- NA
  }
  
  ### compute inflexion point using fitting line
  inflex_point_fit <- sapply(2:(length(x)-1),
                             function(z){
                               slope_1 <- (x[z] - x[1])/(z - 1)
                               intercept_1 <- x[z] - z*slope_1
                               slope_2 <- (x[length(x)] - x[z])/(length(x) - z)
                               intercept_2 <- x[z] - z*slope_2
                               
                               RSS_1 = sum(((1:z * slope_1 + intercept_1) - x[1:z])**2)
                               RSS_2 = sum(((z:length(x) * slope_2 + intercept_2) - x[z:length(x)])**2)
                               
                               res <- RSS_1 + RSS_2
                               
                               return(res)
                             })
  
  inflex_point_fit <- which.min(inflex_point_fit) + 1
  
  ### refilter
  if(is.na(inflex_point_exac)){
    inflex_point <- inflex_point_fit
  }
  else if(inflex_point_fit == inflex_point_exac){
    inflex_point <- inflex_point_fit
  }
  else{
    inflex_point <- c(inflex_point_fit, inflex_point_exac)
    inflex_point <- inflex_point[order(inflex_point)]
    inflex_point <- inflex_point[1]:inflex_point[2]
    
    Ninflex <- abs(inflex_point_fit - inflex_point_exac) + 1
    Npep <- length(x)
    
    if((Ninflex/Npep > 0.75 & Ninflex > 5) | Ninflex > 16){ # 75% of the peptides are considered potential cleavage site --> could have two cleavage around the two points
      # to reduce number of possibilities, checking which method will give greatest mean difference and take its 25% neighbours
      mp_1 <- mean(mp[1:(inflex_point[1]-1)])
      mp_2 <- mean(mp[(inflex_point[Ninflex]+1):Npep])
      Ninflex_25 <- ceiling(Ninflex/4)
      if(mp_1 >= 0 & mp_2 >= 0){
        if(mp_1 < mp_2)
          inflex_point <- inflex_point[1:Ninflex_25]
        else
          inflex_point <- inflex_point[(Ninflex - Ninflex_25 + 1):Ninflex]
      }
      else if(mp_1 < 0 & mp_2 < 0){
        if(mp_1 > mp_2)
          inflex_point <- inflex_point[1:Ninflex_25]
        else
          inflex_point <- inflex_point[(Ninflex - Ninflex_25 + 1):Ninflex]
      }
      else{
        inflex_point <- c(inflex_point[1:Ninflex_25],
                          inflex_point[(Ninflex - Ninflex_25 + 1):Ninflex])
      }
    }
  }
  return(inflex_point)
}

inflex_old <- function(mp){
  ### compute cumulative sum
  if(any(is.na(mp)))
    mp[which(is.na(mp))] <- 0
  
  x <- cumsum(mp)
  
  ### compute 'exact' theoretical inflexion point using second derivative
  x_df1 <- diff(x) # get first 'first derivative'
  # compute 'tendency' from 'first derivative', i.e. if increasing/decreasing
  df <- data.frame(x = 1:length(x_df1),
                   y = x_df1)
  res <- lm(y ~ x, data = df)
  res <- res$coeff[2]
  res <- sign(res)
  
  x_df2 <- diff(x_df1) # get second 'first derivative'
  if(res == 1){
    inflex_point_exac <- which.max(x_df2) + 2 # double diff, so add 2 to index point
  }
  else{
    inflex_point_exac <- which.min(x_df2) + 2 # double diff, so add 2 to index point
  }
  
  if(inflex_point_exac - 2 == length(x_df2)){  # if inflex_point is last one, can't conclude that is cleaved
    inflex_point_exac <- NA
  }
  
  return(inflex_point_exac)
}

### function to sum number of values per peptide N and peptide C
nb_pepval <- function(x){
  if(length(x) < 2)
    return(x)
  
  x <- strsplit(x, "_|\\||;")
  
  y <- x[[1]]
  x <- Reduce(rbind,
              sapply(x,
                     function(y){
                       t(as.data.frame(y))
                     },
                     simplify = FALSE)
  )
  
  x[1,grep("^\\d{1,3}$", y)] <- apply(x[,grep("^\\d{1,3}$", y)], 2, function(z) sum(as.numeric(z)))
  x <- paste(x[1,], collapse = "_")
  x <- gsub("(?<=_\\d{1,3})_(?=\\d{2}C-)", "|", x, perl = TRUE)
  x <- gsub("(?<=_\\d{1,3})_(?=B\\d{1,2}_)", ";", x, perl = TRUE)
  
  return(x)
}

### function to find p-value cutoff
find_cutoff <- function(x,y){
  id <- order(x)
  x <- x[id]
  y <- y[id]
  
  x <- x[which(x < y)]
  x <- ifelse(length(x), x[length(x)], NA)
  return(x)
}

### function to compute cutoff curve
curve <- function(x, cut_neg, cut_pos, cut_p, curvature = 0.1){
  y <- rep(NA, length(x))
  neg <- which(x <= cut_neg)
  pos <- which(x >= cut_pos)
  y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
  y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)
  
  return(y)
}

#### rewriting imprints_cleaved_peptides to allow modification of LMM
# i.e. with or without weights; with or without Rep variable
imprints_cleaved_peptides_LMMeval <- function(data, data_diff = NULL, control = "Vehicle",
                                              min_ValidValue = 0.4, min_peptide = 4,
                                              FDR = 0.01, comb_pv_method = c("george", "fisher", "edgington"),
                                              RESP_score = 0.5, fixed_score_cutoff = FALSE,
                                              categorize = TRUE,
                                              curvature = 0.05, folder_name = "",
                                              use_weights = TRUE, use_RepVar = TRUE # LMM evaluation params
                                              ){
  if(!any(grepl(paste0("_", control, "$"), colnames(data)))){
    message(paste("Error:", control, "is wasn't found in your data !"))
    return()
  }
  
  comb_pv_method <- tolower(comb_pv_method)
  comb_pv_method <- match.arg(comb_pv_method)
  
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_RESP_analysis",
                   ifelse(nchar(folder_name) > 0, "_", ""), folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }
  
  require(limma)
  
  # if any, removing QP columns
  if(any(grepl("^36C_", colnames(data)))){
    data <- data[,-grep("^36C_", colnames(data))]
  }
  
  if(is.null(data_diff)){
    message("No fold-change dataset specified, computing...")
    data_diff <- imprints_sequence_peptides(data, control = control,
                                            dataset_name = deparse(substitute(data))
    )
  }
  
  if(any(grepl("^36C_", colnames(data_diff)))){
    data_diff <- data_diff[,-grep("^36C_", colnames(data_diff))]
  }
  
  if("countNum" %in% colnames(data)){
    data$countNum <- NULL
  }
  if("countNum" %in% colnames(data_diff)){
    data_diff$countNum <- NULL
  }
  
  data_cleaved <- data_diff
  if(any(grepl(paste0("_", control, "$"), colnames(data_cleaved)))){
    data_cleaved <- data_cleaved[,-grep(paste0("_", control, "$"), colnames(data_cleaved))]
  }
  
  # data_cleaved is output from 'imprints_sequence_peptides' --> remove unnecessary informations
  data_cleaved$Master.Protein.Accessions <- data_cleaved$Positions.in.Master.Proteins
  data_cleaved$Positions.in.Master.Proteins <- NULL
  data_cleaved$Annotated.Sequence <- NULL
  data_cleaved$Modifications <- NULL
  colnames(data_cleaved)[1] <- "id" # easier handle
  
  message("Averaging value, summing profiles...")
  # get summed value from all temperature for every peptide + # keeping track of number of non-missing values
  data_cleaved <- data_cleaved %>%
    tidyr::gather("treatment", "value", -id, -description) %>%
    tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(id, description, temp, treatment, biorep) %>% # if peptide modified, mean values from same peptide sequence
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    dplyr::mutate(N = paste(biorep, as.numeric(!is.na(value)), sep = "_")) %>%  # keeping track of number of non-missing values
    dplyr::ungroup() %>% dplyr::group_by(id, description, temp, treatment) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                     Nvalue = paste(N, collapse = ";")) %>%  # keeping track of number of non-missing values
    dplyr::mutate(Nvalue = paste(temp, Nvalue, sep = "-")) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(na.omit(mean_value))/length(mean_value) >= min_ValidValue) %>%  # keeping peptides with more than 40% of valid values
    dplyr::reframe(max_profile = mean_value[which.max(abs(mean_value))],  # get sum value for each peptides --> sum all temperatures values
                   Nvalue = paste(Nvalue, collapse = "|"))   # keeping track of number of non-missing values
  
  message("Filtering and ordering peptides...")
  # separate id in protein and sequence (need to take into account protein groups)
  data_cleaved$sequence <- gsub("(?<= |^)(.{5,6}|A0.{6,8})(?= \\[)", "", data_cleaved$id, perl = TRUE)
  data_cleaved$sequence <- gsub(" ", "", data_cleaved$sequence)
  data_cleaved$sequence <- gsub(";", "; ", data_cleaved$sequence)
  
  data_cleaved$id <- gsub("(?<=\\[)\\d{1,4}-\\d{1,4}(?=\\])", "", data_cleaved$id, perl = TRUE)
  data_cleaved$id <- gsub(" \\[(\\]|\\];)", "", data_cleaved$id)
  data_cleaved$id <- gsub("(?<!;) ", "; ", data_cleaved$id, perl = TRUE)
  
  data_cleaved <- data_cleaved %>%
    dplyr::mutate(sequence = gsub("\\[|\\]", "", sequence)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(max_profile) >= min_peptide)  # only keeping proteins with more than min_peptide peptides
  
  # order data according proteins and sequence
  # ordering according sequence is capital
  data_cleaved$factor <- data_cleaved$id
  data_cleaved <- data_cleaved[order(data_cleaved$id),]
  ord <- sapply(strsplit(gsub(";.*", "", data_cleaved$sequence), "-|~"), # only taking first protein if proteinn group
                function(x) {sum(as.numeric(x))}
  )
  ord <- data.frame(factor = data_cleaved$factor, idx = 1:nrow(data_cleaved), sum.pos = ord) %>%
    dplyr::group_by(factor) %>%
    dplyr::mutate(ord.sum.pos = order(sum.pos),
                  nb.pep = length(ord.sum.pos),
                  final.order = idx[ord.sum.pos])
  
  data_cleaved <- data_cleaved[ord$final.order,]
  data_cleaved$Master.Protein.Accessions <- data_cleaved$factor
  data_cleaved$factor <- NULL
  
  message("Computing cumulative sum to find potential cleaved sites...")
  # get cleaved site
  # if cumulative abundance concave --> first derivative is decreasing globally --> inflexion point is where second derivative is minimum
  # if cumulative abundance convex --> first derivative is increasing globally --> inflexion point is where second derivative is maximum
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::group_modify(~ {
      cleaved_sites = .x$sequence[inflex(.x$max_profile)]
      df <- sapply(cleaved_sites,
                   function(cleaved_site){
                     df <- .x
                     df$cleaved_site <- cleaved_site;
                     df
                   }, simplify = FALSE)
      df <- as.data.frame(Reduce(rbind, df))
      return(df)
    })# get back potential cleaved sites
  
  
  message("Checking cleaved sites, formating...")
  # set peptide position --> before cleaved site = N-terminal side / after cleaved site = C-terminal side
  # if protein group, only taking first protein (take into account overlapping peptides --> tolerance of 2 AA)
  data_cleaved$position_global <- mapply(function(cl, s){
    cl <- as.numeric(cl)
    s <- as.numeric(s)
    
    ifelse(-2 <= s[1] - cl[1] & s[2] - cl[2] <= 2,
           "cleaved",
           ifelse(s[2] <= cl[1], "N",
                  "C")
    )
  },
  strsplit(gsub(";.*", "", data_cleaved$cleaved_site), "~|-"),
  strsplit(gsub(";.*", "", data_cleaved$sequence), "~|-")
  )
  
  # remove cleavage sites giving same positions and keeping widest cleavage site
  data_cleaved <-  data_cleaved %>%
    ungroup() %>%
    group_by(id, description, treatment) %>%
    group_modify(~ {
      simi_pos_seq <- .x %>%
        group_by(cleaved_site) %>%
        summarise(position_sequence = paste(table(position_global), collapse = "-"))
      
      if(any(duplicated(simi_pos_seq$position_sequence))){
        simi_pos_seq <- simi_pos_seq[which(!is.na(
          match(simi_pos_seq$position_sequence,
                simi_pos_seq$position_sequence[which(duplicated(simi_pos_seq$position_sequence))]
          )
        )),
        ] %>%
          ungroup() %>% group_by(position_sequence) %>%
          mutate(size_cleavage = sapply(strsplit(gsub(";.*", "", cleaved_site), "~|-"),
                                        function(s) diff(as.numeric(s))
          ),
          cleavage_torm = size_cleavage != max(size_cleavage)
          )
        
        if(all(!simi_pos_seq$cleavage_torm)){ # meaning none are different from minimum i.e. same size
          # so take lower bound in left and higher bound in right
          simi_pos_seq <- simi_pos_seq %>%
            group_by(position_sequence) %>%
            mutate(cleavage1 = sapply(strsplit(gsub(";.*", "", cleaved_site), "~|-"),
                                      function(s) as.numeric(s)[1]
            ),
            cleavage2 = sapply(strsplit(gsub(";.*", "", cleaved_site), "~|-"),
                               function(s) as.numeric(s)[2]
            )) %>%
            mutate(new_cleaved = ifelse(cleaved_site == cleaved_site[1],
                                        paste0(min(cleavage1), "~", max(cleavage2)),
                                        "torm")
            )
          
          .x <- .x[-which(!is.na(match(.x$cleaved_site,
                                       simi_pos_seq$cleaved_site[which(simi_pos_seq$new_cleaved == "torm")])
          )
          ),]
          
          for(cl in which(simi_pos_seq$new_cleaved != "torm")){
            .x$cleaved_site[which(.x$cleaved_site == simi_pos_seq$cleaved_site[cl])] <- simi_pos_seq$new_cleaved[cl]
          }
          
          return(.x)
        }
        else{ # there exist a  widest cleavage site --> only keeping this one
          simi_pos_seq <- simi_pos_seq$cleaved_site[which(simi_pos_seq$cleavage_torm)]
          .x <- .x[-which(!is.na(match(.x$cleaved_site, simi_pos_seq))),]
          
          return(.x)
        }
      }
      else
        return(.x)
    })
  
  # checking cleaved site, if not assign right position
  data_cleaved <- data_cleaved %>% dplyr::group_by(id, description, treatment, cleaved_site) %>%
    filter(all(c("N", "cleaved", "C") %in% position_global)) %>% # remove proteins not having two sides
    dplyr::mutate(cleaved_site_update = cleaved_site) %>%
    dplyr::group_modify(~ {
      val_N <- .x$max_profile[which(.x$position_global == "N")]
      val_C <- .x$max_profile[which(.x$position_global == "C")]
      cle <- .x$max_profile[which(.x$position_global == "cleaved")]
      non_cleaved <- .x$cleaved_site_update[1]
      npep_cleaved <- length(cle) # in case of overlapping peptides, can have more than one peptide labelled as 'cleaved'
      
      bef_cleaved <- .x$sequence[which(.x$position_global == "cleaved")[1] - 1]
      aft_cleaved <- .x$sequence[which(.x$position_global == "cleaved")[npep_cleaved] + 1]
      
      bef_cleaved <- sapply(strsplit(bef_cleaved, "; ")[[1]],
                            function(x){
                              x <- strsplit(x, "-|~")[[1]]
                              as.numeric(x[2]) + 1
                            }, USE.NAMES = FALSE)
      aft_cleaved <- sapply(strsplit(aft_cleaved, "; ")[[1]],
                            function(x){
                              x <- strsplit(x, "-|~")[[1]]
                              as.numeric(x[1]) - 1
                            }, USE.NAMES = FALSE)
      
      # proportion augmentation/diminution of sd
      pN <- (sd(c(val_N, cle), na.rm = TRUE)/sd(val_N, na.rm = TRUE))
      pC <- (sd(c(val_C, cle), na.rm = TRUE)/sd(val_C, na.rm = TRUE))
      
      if(any(is.na(c(pN, pC)))){ # if only one peptide in N or C --> sd is NA; only one of them can be NA (at least 4 peptides in total)
        # if only one, cleaved site is assign to this group
        .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(is.na(pN), "N", "C")
        non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                              function(x){
                                x <- strsplit(x, "-|~")[[1]]
                                as.numeric(x)
                              }, simplify = FALSE)
        if(is.na(pN)){
          non_cleaved <- sapply(non_cleaved, "[[", 2)
          if(any(non_cleaved != aft_cleaved, na.rm = TRUE)){  # if position doesn't follow
            non_cleaved[which(non_cleaved != aft_cleaved)] <- non_cleaved[which(non_cleaved != aft_cleaved)] + 1
            if(any(non_cleaved > aft_cleaved, na.rm = TRUE)){ # cleaved position and after have common AA position
              non_cleaved[which(non_cleaved > aft_cleaved)] <- aft_cleaved
              aft_cleaved[which(aft_cleaved == non_cleaved)] <- aft_cleaved[which(aft_cleaved == non_cleaved)] + 1
            }
          }
          non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
        }
        else{
          non_cleaved <- sapply(non_cleaved, "[[", 1)
          if(any(non_cleaved != bef_cleaved, na.rm = TRUE)){  # if position doesn't follow
            non_cleaved[which(non_cleaved != bef_cleaved)] <- non_cleaved[which(non_cleaved != bef_cleaved)] - 1
            if(any(non_cleaved < bef_cleaved, na.rm = TRUE)){ # cleaved position and before have common AA position
              bef_cleaved[which(bef_cleaved > non_cleaved)] <- non_cleaved - 1
            }
          }
          non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
        }
        non_cleaved <- paste(non_cleaved, collapse = "; ")
      }
      else if(all(cle >= 0)){ # to be cleaved, peptide should be less abundant, i.e. negative fold changes
        .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(pN <= pC, "N", "C")
        non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                              function(x){
                                x <- strsplit(x, "-|~")[[1]]
                                as.numeric(x)
                              }, simplify = FALSE)
        if(pN <= pC){
          non_cleaved <- sapply(non_cleaved, "[[", 2)
          if(any(non_cleaved != aft_cleaved, na.rm = TRUE)){  # if position doesn't follow
            non_cleaved[which(non_cleaved != aft_cleaved)] <- non_cleaved[which(non_cleaved != aft_cleaved)] + 1
            if(any(non_cleaved > aft_cleaved, na.rm = TRUE)){ # cleaved position and after have common AA position
              non_cleaved[which(non_cleaved > aft_cleaved)] <- aft_cleaved
              aft_cleaved[which(aft_cleaved == non_cleaved)] <- aft_cleaved[which(aft_cleaved == non_cleaved)] + 1
            }
          }
          non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
        }
        else{
          non_cleaved <- sapply(non_cleaved, "[[", 1)
          if(any(non_cleaved != bef_cleaved, na.rm = TRUE)){  # if position doesn't follow
            non_cleaved[which(non_cleaved != bef_cleaved)] <- non_cleaved[which(non_cleaved != bef_cleaved)] - 1
            if(any(non_cleaved < bef_cleaved, na.rm = TRUE)){ # cleaved position and before have common AA position
              bef_cleaved[which(bef_cleaved > non_cleaved)] <- non_cleaved - 1
            }
          }
          non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
        }
        non_cleaved <- paste(non_cleaved, collapse = "; ")
      }
      else{
        if(any(c(pN, pC) < 1)){ # arbitrary, sd augmentation should be more than 100%
          .x$position_global[which(.x$position_global == "cleaved")] <- ifelse(pN <= pC, "N", "C")
          non_cleaved <- sapply(strsplit(non_cleaved, "; ")[[1]],
                                function(x){
                                  x <- strsplit(x, "-|~")[[1]]
                                  as.numeric(x)
                                }, simplify = FALSE)
          
          if(pN <= pC){
            non_cleaved <- sapply(non_cleaved, "[[", 2)
            if(any(non_cleaved != aft_cleaved, na.rm = TRUE)){  # if position doesn't follow
              non_cleaved[which(non_cleaved != aft_cleaved)] <- non_cleaved[which(non_cleaved != aft_cleaved)] + 1
              if(any(non_cleaved > aft_cleaved, na.rm = TRUE)){ # cleaved position and after have common AA position
                non_cleaved[which(non_cleaved > aft_cleaved)] <- aft_cleaved
                aft_cleaved[which(aft_cleaved == non_cleaved)] <- aft_cleaved[which(aft_cleaved == non_cleaved)] + 1
              }
            }
            non_cleaved <- paste(non_cleaved, aft_cleaved, sep = "~")
          }
          else{
            non_cleaved <- sapply(non_cleaved, "[[", 1)
            if(any(non_cleaved != bef_cleaved, na.rm = TRUE)){  # if position doesn't follow
              non_cleaved[which(non_cleaved != bef_cleaved)] <- non_cleaved[which(non_cleaved != bef_cleaved)] - 1
              if(any(non_cleaved < bef_cleaved, na.rm = TRUE)){ # cleaved position and before have common AA position
                bef_cleaved[which(bef_cleaved > non_cleaved)] <- non_cleaved - 1
              }
            }
            non_cleaved <- paste(bef_cleaved, non_cleaved, sep = "~")
          }
          non_cleaved <- paste(non_cleaved, collapse = "; ")
        }
      }
      .x$cleaved_site_update <- non_cleaved
      
      return(.x)
    }) %>%
    dplyr::filter(position_global != "cleaved") %>%
    dplyr::mutate(cleaved_site = cleaved_site_update) %>%
    dplyr::select(-cleaved_site_update) %>% unique()
  
  ### reducing cleavage site candidate
  # remove eventual cleavage site put to last or first position due to many overlapping peptides (rare case)
  # filtering above would then give an NA -->
  cle_torm <- grep("NA", gsub(";.*", "", data_cleaved$cleaved_site))
  if(length(cle_torm))
    data_cleaved <- data_cleaved[-cle_torm,]
  
  # if a protein has 2 or 3 potential cleavage sites and that some are measured peptides,
  # meaning they were decreasing in level and considered as measured cleavage site;
  # then we only keep these as this mean that the two methods used to obtain the inflection point
  # overlapped well and that the measured peptides is very likely to be the actual cleavage site
  data_cleaved <- data_cleaved %>% ungroup() %>%
    dplyr::select(-max_profile, -Master.Protein.Accessions) %>%
    dplyr::group_by(id, description, treatment)  %>%
    dplyr::filter(length(unique(cleaved_site)) == 1 | length(unique(cleaved_site)) > 3 |
                    (length(unique(cleaved_site)) == 2 | length(unique(cleaved_site)) == 3) & (grepl("-", cleaved_site) | all(grepl("~", cleaved_site))))
  
  ## continue to reduce candidates --> remove and summarize cleavage sites creating the same comparison
  data_cleaved <- data_cleaved %>%
    ungroup() %>%
    group_by(id, description, treatment) %>%
    group_modify(~ {
      simi_pos_seq <- .x %>%
        group_by(cleaved_site) %>%
        summarise(position_sequence = paste(table(position_global), collapse = "-"))
      
      if(any(duplicated(simi_pos_seq$position_sequence))){
        simi_pos_seq <- simi_pos_seq[which(!is.na(
          match(simi_pos_seq$position_sequence,
                simi_pos_seq$position_sequence[which(duplicated(simi_pos_seq$position_sequence))]
          )
        )),
        ] %>%
          ungroup() %>% group_by(position_sequence) %>%
          mutate(cleavage1 = sapply(strsplit(gsub(";.*", "", cleaved_site), "~|-"),
                                    function(s) as.numeric(s)[1]
          ),
          cleavage2 = sapply(strsplit(gsub(";.*", "", cleaved_site), "~|-"),
                             function(s) as.numeric(s)[2]
          ),
          are_similar = min(cleavage1 - max(cleavage1)) >= -2 &
            min(cleavage2 - max(cleavage2)) >= -2
          )
        
        if(all(simi_pos_seq$are_similar)){
          simi_pos_seq <- simi_pos_seq %>%
            ungroup() %>% group_by(position_sequence) %>%
            mutate(new_cleaved = ifelse(cleaved_site == cleaved_site[1],
                                        paste0(min(cleavage1), "~", max(cleavage2)),
                                        "torm")
            )
          
          .x <- .x[-which(!is.na(match(.x$cleaved_site,
                                       simi_pos_seq$cleaved_site[which(simi_pos_seq$new_cleaved == "torm")])
          )
          ),]
          
          for(cl in which(simi_pos_seq$new_cleaved != "torm")){
            .x$cleaved_site[which(.x$cleaved_site == simi_pos_seq$cleaved_site[cl])] <- simi_pos_seq$new_cleaved[cl]
          }
          
          return(.x)
        }
        else{
          return(.x)
        }
      }
      else
        return(.x)
    })
  
  ### sum number of values per peptide N and peptide C
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment, cleaved_site, position_global) %>%
    dplyr::summarise(Nvalue = nb_pepval(Nvalue),
                     Npep = length(sequence)
    ) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::mutate(id_cleaved = paste0(id, "_", as.numeric(factor(cleaved_site))))
  
  message("Summing peptides profiles from N-terminal side and C-terminal side from each treatment\n")
  treat <- unique(data_cleaved$treatment)
  new_data_diff <- list()
  res <- list()
  for(t in treat){
    treat_data <- unique(data_cleaved[data_cleaved$treatment == t, c("id_cleaved", "cleaved_site")])
    treat_data_diff <- data
    treat_torm <- treat[!(treat %in% t)]
    if(length(treat_torm)){
      treat_torm <- paste0("_", treat_torm, "$", collapse = "|")
      treat_data_diff <- treat_data_diff[,-grep(treat_torm, colnames(treat_data_diff))]
    }
    
    ## taking into accounts for each potential cleaved sites
    colnames(treat_data_diff)[1] <- "id"
    treat_data_diff <- right_join(treat_data_diff, unique(data_cleaved[data_cleaved$treatment == t,
                                                                       c("id", "id_cleaved")]),
                                  by = "id", relationship = "many-to-many")
    treat_data_diff$id <- treat_data_diff$id_cleaved
    treat_data_diff$id_cleaved <- NULL
    colnames(treat_data_diff)[1] <- "Master.Protein.Accessions"
    ## ids are now updated
    
    directory_toremove <- list.files()
    treat_data_diff <- imprints_sequence_peptides(treat_data_diff,
                                                  proteins = treat_data$id_cleaved,
                                                  sequence = sub("~", "-", treat_data$cleaved_site),
                                                  control = control, dataset_name = t
    )
    directory_toremove <- list.files()[!(list.files() %in% directory_toremove)]
    directory_toremove <- grep(paste0(t, "\\.txt$"), directory_toremove, value = TRUE)
    if(length(directory_toremove) == 1){
      unlink(directory_toremove, recursive = TRUE)
    }
    
    res[[t]] <- treat_data
    ### removing peptides that could be cleaved site
    treat_data_diff <- imprints_remove_peptides(treat_data_diff,
                                                treat_data$id_cleaved,
                                                treat_data$cleaved_site
    )
    # if only one peptide, id on position becomes different than in master, need to modify:
    prid_pos <- gsub(" \\[(\\]|\\])", "",
                     gsub("(?<=\\[)\\d{1,4}(-|~)\\d{1,4}(?=\\])", "",
                          treat_data_diff$Positions.in.Master.Proteins, perl = TRUE)
    )
    if(any(treat_data_diff$Master.Protein.Accessions != prid_pos)){
      prid_new <- treat_data_diff$Master.Protein.Accessions[which(treat_data_diff$Master.Protein.Accessions != prid_pos)]
      prid_new <- gsub(".*; ", "", prid_new)
      prid_pos_tochange <- treat_data_diff$Positions.in.Master.Proteins[which(treat_data_diff$Master.Protein.Accessions != prid_pos)]
      
      treat_data_diff$Positions.in.Master.Proteins[which(treat_data_diff$Master.Protein.Accessions != prid_pos)] <- unlist(
        mapply(function(x, n){
          if(grepl("; ", x)){ # protein group and/or different possible peptides
            if(grepl("; \\[", x)){ # different possible peptides
              x <- gsub(sub("_.*", "", n),
                        n, x)
            }
            else{ # protein group
              x <- gsub(paste0("; ", sub("_.*", "", n), ".*"),
                        "", x)
              x <- paste0(x, "; ", n)
            }
          }
          
          else{
            x <- gsub(sub("_.*", "", n),
                      n, x)
          };
          x
        },
        prid_pos_tochange, prid_new,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
      )
    }
    
    treat_data_diff <- treat_data_diff[,-grep(paste0("_", control), colnames(treat_data_diff))]
    treat_data_diff$Master.Protein.Accessions <- treat_data_diff$Positions.in.Master.Proteins
    treat_data_diff$Modifications <- NULL
    treat_data_diff$Annotated.Sequence <- NULL
    treat_data_diff$Positions.in.Master.Proteins <- NULL
    colnames(treat_data_diff)[1] <- "id"
    
    treat_data_diff <- treat_data_diff %>%
      tidyr::gather("key", "value", -id, -description) %>%
      tidyr::separate(key, into = c("temperature", "biorep", "treatment"), sep = "_") %>%
      dplyr::mutate(protein = gsub("(?<!;) ", "; ", gsub(" \\[(\\]|\\];)", "", gsub("(?<=\\[)\\d{1,4}(-|~)\\d{1,4}(?=\\])", "", id, perl = TRUE)), perl = TRUE),
                    position = gsub(";", "; ", gsub(" ", "", gsub("(?<=; |^)(.{5,6}|A0.{6,8})(|_\\d{1,3})(?= \\[)", "", id, perl = TRUE)))
      ) %>%
      dplyr::select(-id)
    treat_data_diff <- treat_data_diff[,c("protein", "position", "description", "temperature", "biorep", "treatment", "value")]
    new_data_diff[[t]] <- treat_data_diff
  }
  new_data_diff <- as.data.frame(Reduce(rbind, new_data_diff))
  # updating ids
  data_cleaved$id <- data_cleaved$id_cleaved
  data_cleaved$id_cleaved <- NULL
  
  ##### computing p-values for each protein --> H0 = protein is not cleaved
  message("\nReshaping data and computing experimental design")
  new_data_diff$position_global <- unlist(lapply(strsplit(gsub("\\[|\\]", "", gsub(";.*", "", new_data_diff$position)),
                                                          "~|-"),
                                                 function(s) sum(as.numeric(s)))
  )
  new_data_diff <- new_data_diff %>%
    dplyr::group_by(protein, description, treatment) %>%
    dplyr::mutate(position_global = ifelse(position_global == min(position_global),
                                           "N", "C")) %>%
    dplyr::select(-position)
  
  new_data_diff$temperature <- factor(new_data_diff$temperature)
  new_data_diff$biorep <- factor(new_data_diff$biorep)
  new_data_diff$position_global <- factor(new_data_diff$position_global,
                                          levels = c("N", "C"))
  
  ntemp <- nlevels(new_data_diff$temperature)
  nrep <- nlevels(new_data_diff$biorep)
  # npos = 2 ; N or C
  
  ### limma design
  # design
  position <- rep(levels(new_data_diff$position_global), ntemp*nrep)
  temperature <- rep(levels(new_data_diff$temperature), each=nrep*2)
  replicate <- rep(rep(levels(new_data_diff$biorep), each = 2), ntemp)
  subj <- paste0(position, replicate)
  # g features points of interest to build contrasts
  g <- factor(paste0(position, temperature),
              levels = c(paste0(levels(new_data_diff$position_global)[1],
                                levels(new_data_diff$temperature)
              ),
              paste0(levels(new_data_diff$position_global)[2],
                     levels(new_data_diff$temperature)
              )
              )
  )
  
  if(use_RepVar)
    design <- model.matrix(~ 0 + g + replicate)
  else
    design <- model.matrix(~ 0 + g)
  
  colnames(design)[1:nlevels(g)] <- levels(g)
  # reshape data and order columns
  new_data_diff_matrix <- new_data_diff %>%
    tidyr::unite("key", position_global, biorep, temperature, sep = "_") %>%
    tidyr::spread(key, value) %>%
    dplyr::mutate(description = unname(sapply(description, IMPRINTS.CETSA.app:::getGeneName))) %>%
    tidyr::unite("id", protein, description, sep = "_")
  new_data_diff_matrix <- new_data_diff_matrix[,c("id", "treatment", paste(position, replicate, temperature, sep = "_"))]
  
  ##### Get weight matrix for lmFit
  if(use_weights){
    message("Computing weight matrix from number of measurements")
    weights <- data_cleaved
    weights$description <- unname(sapply(weights$description, IMPRINTS.CETSA.app:::getGeneName))
    weights$cleaved_site <- NULL
    weights$Npep <- NULL
    
    ## extract number of measurements
    weights <- weights %>%
      dplyr::group_by(id, description, treatment, position_global) %>%
      dplyr::group_modify(~ {
        x <- lapply(strsplit(.x$Nvalue, "\\|")[[1]], function(n) strsplit(n, "-")[[1]])
        
        temp <- unlist(lapply(x, "[[", 1))
        biorep <- strsplit(unlist(lapply(x, "[[", 2)), ";")
        
        x <- mapply(function(b, t){
          b <- strsplit(b, "_")
          b <- t(as.data.frame(b))
          rownames(b) <- 1:nrow(b)
          colnames(b) <- c("biorep", "Nv")
          b <- as.data.frame(b)
          b$Nv <- as.numeric(b$Nv)
          b$temperature <- t;
          b
        },
        biorep, temp, SIMPLIFY = FALSE)
        
        x <- as.data.frame(Reduce(rbind, x))
        return(x)
      })
    
    weights <- weights %>%
      tidyr::unite("id", id, description, sep = "_") %>%
      tidyr::unite("key", position_global, biorep, temperature, sep = "_") %>%
      tidyr::spread(key, Nv)
    weights <- weights[,colnames(new_data_diff_matrix)]
  }
  
  message("Computing p-values for each temperatures for each treatment...")
  res <- list()
  for(tr in treat){
    X <- new_data_diff_matrix[new_data_diff_matrix$treatment == tr,] %>%
      dplyr::ungroup() %>% dplyr::select(-treatment) %>%
      tibble::column_to_rownames("id") %>%
      as.matrix()
    
    if(use_weights){
      Xw <- weights[weights$treatment == tr,] %>%
        dplyr::ungroup() %>% dplyr::select(-treatment) %>%
        tibble::column_to_rownames("id") %>%
        as.matrix()
    }
    else{
      Xw <- NULL
    }
    
    ### Add blocking on subject in duplicateCorrelation().
    corfit <- duplicateCorrelation(X, design, block=subj)
    ### Fit our model with weigths (not normalized)
    fitW <- lmFit(X, design, weights = Xw,
                  block=subj, correlation=corfit$consensus)
    
    temperature <- levels(new_data_diff$temperature)
    res[[tr]] <- list()
    for(tmp in temperature){
      # How do the two profiles differ in their response over temperature ?
      con.1 <- makeContrasts(contrasts = paste0("C", tmp, " - N", tmp),
                             levels = design)
      fitcon.1 <- contrasts.fit(fitW, con.1)
      fitcon.1 <- eBayes(fitcon.1)
      res[[tr]][[tmp]] <- topTable(fitcon.1, number=nrow(X), adjust.method="BH")
    }
    
    res[[tr]] <- as.data.frame(Reduce(rbind,
                                      mapply(function(x, n){
                                        x$temperature <- n
                                        x <- x %>%
                                          tibble::rownames_to_column("id") %>%
                                          tidyr::separate(id, into = c("id", "Gene"), sep = "(?<=_\\d{1,3})_");
                                        x
                                      },
                                      res[[tr]], names(res[[tr]]),
                                      SIMPLIFY = FALSE)
    )
    )
    
    res[[tr]]$treatment <- tr
  }
  res <- as.data.frame(Reduce(rbind, res))
  
  message("Finding most likely proteins with RESP effect")
  res <- res %>%
    dplyr::mutate(adj.P.Val = tidyr::replace_na(adj.P.Val, 0.999999))
  
  if(comb_pv_method == "george"){
    res <- res %>%
      dplyr::group_by(id, Gene, treatment) %>%
      dplyr::summarise(combined_pvalue = ifelse(length(na.omit(logFC)),
                                                ifelse(ntemp > 1, metap::logitp(adj.P.Val)$p, adj.P.Val), # George's method
                                                NA),
                       maxFC = ifelse(length(na.omit(logFC)),
                                      temperature[which(abs(logFC) == max(abs(logFC), na.rm = T))],
                                      NA)
      )
  }
  else if(comb_pv_method == "fisher"){
    res <- res %>%
      dplyr::group_by(id, Gene, treatment) %>%
      dplyr::summarise(combined_pvalue = ifelse(length(na.omit(logFC)),
                                                ifelse(ntemp > 1, metap::sumlog(adj.P.Val)$p, adj.P.Val), # Fisher's method
                                                NA),
                       maxFC = ifelse(length(na.omit(logFC)),
                                      temperature[which(abs(logFC) == max(abs(logFC), na.rm = T))],
                                      NA)
      )
  }
  else if(comb_pv_method == "edgington"){
    res <- res %>%
      dplyr::group_by(id, Gene, treatment) %>%
      dplyr::summarise(combined_pvalue = ifelse(length(na.omit(logFC)),
                                                ifelse(ntemp > 1, metap::sump(adj.P.Val)$p, adj.P.Val), # Edgington’s method
                                                NA),
                       maxFC = ifelse(length(na.omit(logFC)),
                                      temperature[which(abs(logFC) == max(abs(logFC), na.rm = T))],
                                      NA)
      )
  }
  
  # compute RESP-score
  res <- res %>%
    dplyr::group_by(id, Gene, treatment) %>%
    dplyr::group_modify(~ {
      if(!is.na(.x$maxFC)){
        # getting FC
        score <- new_data_diff_matrix[which(new_data_diff_matrix$id == paste0(.x$id, "_", .x$Gene) & new_data_diff_matrix$treatment == .x$treatment),
                                      grep(paste0("_", .x$maxFC, "$"), colnames(new_data_diff_matrix))]
        
        # computing score
        score <- c(mean(as.numeric(score[1,grep("^N_", colnames(score))]), na.rm = TRUE),
                   mean(as.numeric(score[1,grep("^C_", colnames(score))]), na.rm = TRUE))
        sign_score <- unique(sign(score))
        diff_score <- diff(score)
        if(length(unique(sign_score)) == 1){ # same sign
          resp_score <- diff_score/max(c(0.1, min(abs(score))))
          # ratio give the score where 0.1 is a constant that prevent exploding score when min approach 0
        }
        else{ # one is neg and one is pos
          resp_score <- max(abs(score))*diff_score/0.1  # also use the same constant
          # multiplying by the maximum fold change prevent keepin low profile that would have opposite sign
        }
      }
      else{
        resp_score <- NA
      }
      
      score <- .x
      score <- score[,-grep("^id$|^Gene$|^treatment$", colnames(score))]
      score$resp_score <- resp_score
      return(score)
    }, .keep = TRUE)
  
  # selecting cleaved site with smallest p-value and greatest RESP-score
  res <- res %>%
    dplyr::ungroup() %>% dplyr::select(-maxFC) %>%
    filter(!is.na(combined_pvalue)) %>%
    dplyr::mutate(id2 = sub("_.*", "", id)) %>%
    dplyr::group_by(id2, Gene, treatment) %>%
    mutate(rnk = order(order(combined_pvalue)) + order(order(abs(resp_score), decreasing = TRUE)))%>%
    filter(rnk == min(rnk)) %>%
    dplyr::summarise(id = id[which.min(combined_pvalue)], # if more than one cleavage site creates the same rank, take the one with smallest p-value
                     combined_pvalue = min(combined_pvalue, na.rm = TRUE),
                     resp_score = resp_score[which.min(combined_pvalue)]) %>%
    dplyr::ungroup() %>% dplyr::select(-id2)
  
  # z-core normalize RESP-score
  res <- res %>% dplyr::group_by(treatment) %>%
    dplyr::mutate(resp_score = (resp_score - mean(resp_score, na.rm = TRUE))/sd(resp_score, na.rm = TRUE))
  
  # compute cutoffs
  if(fixed_score_cutoff){
    cutoff <- res %>% dplyr::ungroup() %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(combined_pvalue))/length(combined_pvalue))*FDR) %>%
      dplyr::summarise(pval = find_cutoff(combined_pvalue, BH),
                       resp_score_pos = RESP_score,
                       resp_score_neg = -RESP_score)
  }
  else{
    cutoff <- res %>% dplyr::ungroup() %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(combined_pvalue))/length(combined_pvalue))*FDR) %>%
      dplyr::summarise(pval = find_cutoff(combined_pvalue, BH),
                       resp_score_pos = RESP_score + median(abs(resp_score)[which(combined_pvalue < quantile(combined_pvalue, 0.5, na.rm = TRUE))], na.rm = TRUE),
                       resp_score_neg = -RESP_score - median(abs(resp_score)[which(combined_pvalue < quantile(combined_pvalue, 0.5, na.rm = TRUE))], na.rm = TRUE))
  }
  
  # saving obtained cutoffs
  cutoff_file <- paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "cutoff.txt")
  readr::write_tsv(cutoff, cutoff_file)
  extra_info <- paste0("\nParameters: \nRESP z-score cutoff=", RESP_score, ", FDR=", FDR*100, "%, method to combine p-values=", comb_pv_method, ", curvature=", curvature,
                       " \nMinimum number of peptides=", min_peptide, ", Minimum proportion of non-missing values=", min_ValidValue)
  write(extra_info, cutoff_file, sep = "\n", append = TRUE)
  
  res <- res %>% dplyr::group_by(id, Gene, treatment) %>%
    dplyr::mutate(curve = curve(resp_score,
                                cutoff$resp_score_neg[which(cutoff$treatment == treatment)],
                                cutoff$resp_score_pos[which(cutoff$treatment == treatment)],
                                cutoff$pval[which(cutoff$treatment == treatment)],
                                curvature = curvature),
                  criteria_curve = -log10(combined_pvalue) >= curve
    )
  res$criteria_curve <- tidyr::replace_na(res$criteria_curve, FALSE)
  
  n_treat <- length(treat)
  df_curve <- data.frame(resp_score = rep(seq(-max(abs(res$resp_score), na.rm = TRUE),
                                              max(abs(res$resp_score), na.rm = TRUE), 0.005),
                                          n_treat)
  )
  df_curve$treatment <- rep(treat, each = nrow(df_curve)/n_treat)
  df_curve <- df_curve %>% dplyr::group_by(treatment, rownames(df_curve)) %>%
    dplyr::mutate(curve = curve(resp_score,
                                cutoff$resp_score_neg[which(cutoff$treatment == treatment)],
                                cutoff$resp_score_pos[which(cutoff$treatment == treatment)],
                                cutoff$pval[which(cutoff$treatment == treatment)],
                                curvature = curvature
    )
    )
  
  message("Creating and saving plot")
  g_h <- ggplot(res, aes(resp_score, -log10(combined_pvalue), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = resp_score, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(res$combined_pvalue), na.rm = TRUE))) +
    xlim(c(-max(abs(res$resp_score), na.rm = TRUE),
           max(abs(res$resp_score), na.rm = TRUE))
    ) +
    labs(title = "RESP volcano plot",
         y = "-log10(combined p-value)",
         x = "RESP z-score") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    facet_wrap(~treatment) +
    ggrepel::geom_label_repel(data = res[res$criteria_curve,],
                              aes(resp_score, -log10(combined_pvalue), label = Gene),
                              show.legend = FALSE, min.segment.length = 0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))
  
  ggsave(paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "RESP_hits_plot.png"),
         plot = g_h, device = "png",
         width = 16, height = 9)
  
  message("Preparing and saving results")
  data_cleaved$description <- unname(sapply(data_cleaved$description, IMPRINTS.CETSA.app:::getProteinName))
  data_cleaved$position_global <- paste0(data_cleaved$position_global, "-term")
  data_cleaved <- data_cleaved %>%
    tidyr::pivot_wider(id_cols = c("id", "description", "treatment", "cleaved_site"),
                       names_from = position_global,
                       values_from = c("Nvalue", "Npep"))
  
  # saving whole results
  resp <- res[,c("id", "Gene", "treatment", "combined_pvalue", "resp_score", "criteria_curve")]
  colnames(resp)[c(ncol(resp)-1, ncol(resp))] <- c("RESP_score", "RESP_hit")
  resp <- dplyr::left_join(resp, data_cleaved, by = c("id", "treatment"),
                           relationship = "one-to-many")
  resp <- resp[order(resp$combined_pvalue),]
  resp$id <- sub("_.*", "", resp$id)
  openxlsx::write.xlsx(resp, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "RESP_analysis_full.xlsx"))
  
  # saving summary
  resp_summary <- res[res$criteria_curve,c("id", "Gene", "treatment", "combined_pvalue", "resp_score")]
  colnames(resp_summary)[ncol(resp_summary)] <- "RESP_score"
  resp_summary <- dplyr::left_join(resp_summary, data_cleaved, by = c("id", "treatment"),
                                   relationship = "one-to-many")
  resp_summary <- resp_summary[order(resp_summary$combined_pvalue),]
  resp_summary$id <- sub("_.*", "", resp_summary$id)
  
  if(categorize & nrow(resp_summary)){
    message("Categorizing...")
    resp_summary <- imprints_categorize_peptides(data, resp_summary, control, save_xlsx = FALSE)
  }
  
  openxlsx::write.xlsx(resp_summary, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "RESP_hits_summary.xlsx"))
  message("Done !")
  return(resp_summary)
}

#### rewriting imprints_cleaved_peptides function to only apply R-squared method as explained in Ramos et al.
imprints_cleaved_peptides_R2method <- function(data_diff, control = "Vehicle",
                                               R2_cutoff = 0.9, min_ValidValue = 0.4, min_peptide = 4
                                               ){
  if(!any(grepl(paste0("_", control, "$"), colnames(data_diff)))){
    message(paste("Error:", control, "is wasn't found in your data !"))
    return()
  }
  
  if(any(grepl("^36C_", colnames(data_diff)))){
    data_diff <- data_diff[,-grep("^36C_", colnames(data_diff))]
  }
  
  if("countNum" %in% colnames(data_diff)){
    data_diff$countNum <- NULL
  }
  
  data_cleaved <- data_diff
  if(any(grepl(paste0("_", control, "$"), colnames(data_cleaved)))){
    data_cleaved <- data_cleaved[,-grep(paste0("_", control, "$"), colnames(data_cleaved))]
  }
  
  # data_cleaved is output from 'imprints_sequence_peptides' --> remove unnecessary informations
  data_cleaved$Master.Protein.Accessions <- data_cleaved$Positions.in.Master.Proteins
  data_cleaved$Positions.in.Master.Proteins <- NULL
  data_cleaved$Annotated.Sequence <- NULL
  data_cleaved$Modifications <- NULL
  colnames(data_cleaved)[1] <- "id" # easier handle
  
  message("Averaging value, summing profiles...")
  # get summed value from all temperature for every peptide + # keeping track of number of non-missing values
  data_cleaved <- data_cleaved %>%
    tidyr::gather("treatment", "value", -id, -description) %>%
    tidyr::separate(treatment, into = c("temp", "biorep", "treatment"), sep = "_") %>%
    dplyr::group_by(id, description, temp, treatment, biorep) %>% # if peptide modified, mean values from same peptide sequence
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    dplyr::mutate(N = paste(biorep, as.numeric(!is.na(value)), sep = "_")) %>%  # keeping track of number of non-missing values
    dplyr::ungroup() %>% dplyr::group_by(id, description, temp, treatment) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                     Nvalue = paste(N, collapse = ";")) %>%  # keeping track of number of non-missing values
    dplyr::mutate(Nvalue = paste(temp, Nvalue, sep = "-")) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(na.omit(mean_value))/length(mean_value) >= min_ValidValue) %>%  # keeping peptides with more than 40% of valid values
    dplyr::reframe(sum_profile = sum(mean_value, na.rm = TRUE),  # get sum value for each peptides (!! not max !!) --> sum all temperatures values
                   Nvalue = paste(Nvalue, collapse = "|"))   # keeping track of number of non-missing values
  
  message("Filtering and ordering peptides...")
  # separate id in protein and sequence (need to take into account protein groups)
  data_cleaved$sequence <- gsub("(?<= |^)(.{5,6}|A0.{6,8})(?= \\[)", "", data_cleaved$id, perl = TRUE)
  data_cleaved$sequence <- gsub(" ", "", data_cleaved$sequence)
  data_cleaved$sequence <- gsub(";", "; ", data_cleaved$sequence)
  
  data_cleaved$id <- gsub("(?<=\\[)\\d{1,4}-\\d{1,4}(?=\\])", "", data_cleaved$id, perl = TRUE)
  data_cleaved$id <- gsub(" \\[(\\]|\\];)", "", data_cleaved$id)
  data_cleaved$id <- gsub("(?<!;) ", "; ", data_cleaved$id, perl = TRUE)
  
  data_cleaved <- data_cleaved %>%
    dplyr::mutate(sequence = gsub("\\[|\\]", "", sequence)) %>%
    dplyr::ungroup() %>% dplyr::group_by(id, description, treatment) %>%
    dplyr::filter(length(sum_profile) >= min_peptide)  # only keeping proteins with more than min_peptide peptides
  
  # order data according proteins and sequence
  # ordering according sequence is capital
  data_cleaved$factor <- data_cleaved$id
  data_cleaved <- data_cleaved[order(data_cleaved$id),]
  ord <- sapply(strsplit(gsub(";.*", "", data_cleaved$sequence), "-|~"), # only taking first protein if protein group
                function(x) {sum(as.numeric(x))}
  )
  ord <- data.frame(factor = data_cleaved$factor, idx = 1:nrow(data_cleaved), sum.pos = ord) %>%
    dplyr::group_by(factor) %>%
    dplyr::mutate(ord.sum.pos = order(sum.pos),
                  nb.pep = length(ord.sum.pos),
                  final.order = idx[ord.sum.pos])
  
  data_cleaved <- data_cleaved[ord$final.order,]
  data_cleaved$Master.Protein.Accessions <- data_cleaved$factor
  data_cleaved$factor <- NULL
  
  message("Computing cumulative sum to find potential cleaved sites...")
  # get cleaved site
  # if cumulative abundance concave --> first derivative is decreasing globally --> inflexion point is where second derivative is minimum
  # if cumulative abundance convex --> first derivative is increasing globally --> inflexion point is where second derivative is maximum
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::mutate(R2 = cor(1:length(sum_profile), cumsum(sum_profile))**2) 
  all_rsquared <<- data_cleaved
  data_cleaved <- data_cleaved %>%
    dplyr::ungroup() %>% dplyr::filter(R2 <= R2_cutoff) %>%
    dplyr::group_by(id, description, treatment) %>%
    dplyr::group_modify(~ {
      cleaved_sites = .x$sequence[inflex_old(.x$sum_profile)] # using old method: only take exact inflexion point
      df <- sapply(cleaved_sites,
                   function(cleaved_site){
                     df <- .x
                     df$cleaved_site <- cleaved_site;
                     df
                   }, simplify = FALSE)
      df <- as.data.frame(Reduce(rbind, df))
      return(df)
    }) %>%
    dplyr::filter(!is.na(cleaved_site))
  
  
  message("Checking cleaved sites, formating...")
  # set peptide position --> before cleaved site = N-terminal side / after cleaved site = C-terminal side
  # if protein group, only taking first protein (take into account overlapping peptides --> tolerance of 2 AA)
  data_cleaved$position_global <- mapply(function(cl, s){
    cl <- as.numeric(cl)
    s <- as.numeric(s)
    
    ifelse(-2 <= s[1] - cl[1] & s[2] - cl[2] <= 2,
           "cleaved",
           ifelse(s[2] <= cl[1], "N",
                  "C")
    )
  },
  strsplit(gsub(";.*", "", data_cleaved$cleaved_site), "~|-"),
  strsplit(gsub(";.*", "", data_cleaved$sequence), "~|-")
  )
  
  # removing peptide being potiential cleavage site
  data_cleaved <- data_cleaved[which(data_cleaved$position_global != "cleaved"),]
  
  # removing proteins with only one side
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment) %>%
    filter(all(c("C", "N") %in% position_global))
  
  # removing unnecessary columns
  data_cleaved$Master.Protein.Accessions <- NULL
  data_cleaved$sum_profile <- NULL
  
  ### sum number of values per peptide N and peptide C
  data_cleaved <- data_cleaved %>% dplyr::ungroup() %>%
    dplyr::group_by(id, description, treatment, cleaved_site, position_global, R2) %>%
    dplyr::summarise(Nvalue = nb_pepval(Nvalue),
                     Npep = length(sequence)
    ) 
  
  message("Preparing and saving results")
  data_cleaved$Gene <- unname(sapply(data_cleaved$description, IMPRINTS.CETSA.app:::getGeneName))
  data_cleaved$description <- unname(sapply(data_cleaved$description, IMPRINTS.CETSA.app:::getProteinName))
  data_cleaved$position_global <- paste0(data_cleaved$position_global, "-term")
  data_cleaved <- data_cleaved %>%
    tidyr::pivot_wider(id_cols = c("id", "description", "Gene", "treatment", "cleaved_site", "R2"),
                       names_from = position_global,
                       values_from = c("Nvalue", "Npep"))
  
  return(data_cleaved)
}


### read data
vene6h_post_norm <- readr::read_tsv("Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax6h_proteins_postnorm.txt")
vene6h_post_norm <- vene6h_post_norm[,-grep("_Aza|_Dec", colnames(vene6h_post_norm))]
vene6h_diff <- IMPRINTS.CETSA::imprints_caldiff_f(vene6h_post_norm, "Vehicle6h")
# moving to Data folder
file.rename(list.files(pattern = "\\d{6}_\\d{4}_vene6h_post_norm_imprints_caldiff", recursive = TRUE),
            "Data/Venetoclax&AT-IAP_Ramosetal/Venetoclax6h_proteins_diff.txt")

## rerun old method using R2
respect_old_vene6h <- imprints_cleaved_peptides_R2method(vene6h_diff, control = "Vehicle6h", R2_cutoff = 0.9,
                                                         min_peptide = 5, min_ValidValue = 0.5)

## rerun RESPect algorithm
LMMeval_vene6h <- imprints_cleaved_peptides_LMMeval(vene6h_post_norm, vene6h_diff, control = "Vehicle6h", folder_name = "venetoclax6h",
                                                    fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                                    FDR = 0.01, comb_pv_method = "george",
                                                    min_peptide = 5, min_ValidValue = 0.5)
file.rename(list.files(pattern = "\\d{6}_\\d{4}_RESP_analysis_venetoclax6h"),
            "Results/RESP_analysis_venetoclax6h")


## comparing results (SFig 4A)
LMMeval_vene6h_full <- openxlsx::read.xlsx(list.files(path = "Results/RESP_analysis_venetoclax6h", 
                                                      pattern = "analysis_tab", full.names = TRUE))
cutoff <- readr::read_tsv(list.files(path = "Results/RESP_analysis_venetoclax6h", 
                                     pattern = "cutoff", full.names = TRUE),
                          n_max = 1)

LMMeval_vene6h_full$treatment <- "Venetoclax 6h"
cutoff$treatment <- "Venetoclax 6h"

LMMeval_vene6h_full$R2_hit <- LMMeval_vene6h_full$Gene %in% respect_old_vene6h$Gene
LMMeval_vene6h_full$R2_RESPect_comparison <- factor(ifelse(LMMeval_vene6h_full$RESP_hit & LMMeval_vene6h_full$R2_hit,
                                                           "RESPect and R-squared",
                                                           ifelse(LMMeval_vene6h_full$RESP_hit, "RESPect",
                                                                  ifelse(LMMeval_vene6h_full$R2_hit, "R-squared",
                                                                         "None")
                                                                  )
                                                           ),
                                                    levels = c("RESPect and R-squared", "RESPect",
                                                               "R-squared", "None")
                                                    )

df_curve <- data.frame(RESP_score = seq(-max(abs(LMMeval_vene6h_full$RESP_score), na.rm = TRUE),
                                        max(abs(LMMeval_vene6h_full$RESP_score), na.rm = TRUE), 0.001)
                       )
df_curve$treatment <- "Venetoclax 6h"
df_curve <- df_curve %>% 
  dplyr::mutate(curve = IMPRINTS.CETSA.app:::curve(RESP_score, 
                                                   cutoff$resp_score_neg, cutoff$resp_score_pos,
                                                   cutoff$pval))

ggplot(LMMeval_vene6h_full, aes(RESP_score, -log10(combined_pvalue),
                                color = R2_RESPect_comparison)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("RESPect and R-squared" = "#2c6e49", "RESPect" = "#007ea7",
                                "R-squared" = "#8d0801", "None" = "grey70"),
                     labels = paste0(levels(LMMeval_vene6h_full$R2_RESPect_comparison), 
                                     " (n=", sapply(levels(LMMeval_vene6h_full$R2_RESPect_comparison), 
                                                    function(g) length(unique(LMMeval_vene6h_full$Gene[which(LMMeval_vene6h_full$R2_RESPect_comparison == g)]))
                                     ),
                                     ")")) +
  geom_line(data = df_curve, aes(x = RESP_score, y = curve), linetype = "dashed", color = "black") +
  ylim(c(0, max(-log10(LMMeval_vene6h_full$combined_pvalue), na.rm = TRUE))) +
  xlim(c(-max(abs(LMMeval_vene6h_full$RESP_score), na.rm = TRUE), max(abs(LMMeval_vene6h_full$RESP_score), na.rm = TRUE))) +
  labs(title = "RESPect Volcano plot: evaluate and compare R-squared baseline with RESPect",
       y = "-log10(Combined p-value)", x = "RESP score") +
  facet_wrap(~treatment) +
  ggrepel::geom_label_repel(data = LMMeval_vene6h_full[which(LMMeval_vene6h_full$RESP_hit |
                                                               (LMMeval_vene6h_full$R2_hit & 
                                                                  (LMMeval_vene6h_full$combined_pvalue <= 0.01 | abs(LMMeval_vene6h_full$RESP_score) >= 4)
                                                                )
                                                             ),],
                            aes(RESP_score, -log10(combined_pvalue), label = Gene), 
                            show.legend = FALSE, min.segment.length = 0, size = 6) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", color = "black", size = 18),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))

ggsave("Figures/RESPectVolcano_R2eval-compare_Venetoclax6h.png", dpi = 300, device = "png", 
       width = 22, height = 14)


### evaluation of the weighted Linear Mixed Model (LMM): FC ~ PoT + Rep  with blocking variable on PoRep and number of measured peptides as weights
LMMeval_vene6h_noWeights <- imprints_cleaved_peptides_LMMeval(vene6h_post_norm, xd, control = "Vehicle6h", folder_name = "venetoclax6h_noWeights",
                                                              fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                                              FDR = 0.01, comb_pv_method = "george",
                                                              min_peptide = 5, min_ValidValue = 0.5,
                                                              use_weights = FALSE)

LMMeval_vene6h_noRep <- imprints_cleaved_peptides_LMMeval(vene6h_post_norm, xd, control = "Vehicle6h", folder_name = "venetoclax6h_noRep",
                                                          fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                                          FDR = 0.01, comb_pv_method = "george",
                                                          min_peptide = 5, min_ValidValue = 0.5,
                                                          use_RepVar = FALSE)

LMMeval_vene6h_noRep_noWeights <- imprints_cleaved_peptides_LMMeval(vene6h_post_norm, xd, control = "Vehicle6h", folder_name = "venetoclax6h_noRep_noWeights",
                                                                    fixed_score_cutoff = TRUE, RESP_score = 0.5, 
                                                                    FDR = 0.01, comb_pv_method = "george",
                                                                    min_peptide = 5, min_ValidValue = 0.5,
                                                                    use_weights = FALSE, use_RepVar = FALSE)

l <- com_protein_loop(list("Normal" = LMMeval_vene6h$id, 
                           "NoWeights" = LMMeval_vene6h_noWeights$id,
                           "NoRep" = LMMeval_vene6h_noRep$id, 
                           "NoRep_NoWeights" = LMMeval_vene6h_noRep_noWeights$id)) 

# |Normal                                       |  2|
# |NoWeights                                    | 14|
# |Normal & NoWeights                           |  7|
# |NoWeights & NoRep_NoWeights                  |  5|
# |Normal & NoWeights & NoRep                   |  1|
# |Normal & NoWeights & NoRep_NoWeights         | 12|
# |Normal & NoWeights & NoRep & NoRep_NoWeights | 70|


### print venn diagram of these genes (SFig 4B)
venn <- VennDiagram::venn.diagram(list("Normal" = LMMeval_vene6h$Gene, 
                                       "NoWeights" = LMMeval_vene6h_noWeights$Gene,
                                       "NoRep" = LMMeval_vene6h_noRep$Gene, 
                                       "NoRep_NoWeights" = LMMeval_vene6h_noRep_noWeights$Gene), 
                                  filename = NULL, disable.logging = TRUE, 
                                  fill = c("#e63946", "#f1faee", "#a8dadc", "#457b9d"), 
                                  lwd = 2, cex = 3)
png("Figures/LMMeval_nHits_VennDiagram.png", 
    width=12, height=12, units = "in", res=300)
grid::grid.newpage()
grid::grid.draw(venn)
dev.off()


## Not found by Normal --> NoWeights and NoRep_NoWeights
LMMeval_summary <-  sapply(list("Normal" = LMMeval_vene6h, 
                                "NoWeights" = LMMeval_vene6h_noWeights,
                                "NoRep" = LMMeval_vene6h_noRep,
                                "NoRep_NoWeights" = LMMeval_vene6h_noRep_noWeights),
                           function(x) c(table(x$category), "Total" = nrow(x)))

# total population:
nrow(LMMeval_vene6h_full) # 1706
# Prevalence:
knitr::kable(round(100*sort(LMMeval_summary["Total",]/nrow(LMMeval_vene6h_full), decreasing = TRUE), 2))
# ~Positive Predictive Value (PPV)
knitr::kable(round(100*sort(LMMeval_summary["RESP",]/(LMMeval_summary["Unconfident",] + LMMeval_summary["RESP",]), decreasing = TRUE), 2))

knitr::kable(round(100*sort(LMMeval_summary["RESP",]/LMMeval_summary["Total",], decreasing = TRUE), 2))
knitr::kable(round(100*sort(LMMeval_summary["Unconfident",]/LMMeval_summary["Total",]), 2))


