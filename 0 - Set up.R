# R version : 4.5.0

# Checking if devtools is installed, if not installing it
if(!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
} 

# Checking if BiocManager is installed, if not installing it
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager") 
}

# install needed bioconductor packages
BiocManager::install(c("limma", "arrayQualityMetrics",
                       "STRINGdb", "clusterProfiler", 
                       "enrichplot", "org.Hs.eg.db", "GOSemSim",
                       "multtest", "vsn"))

# Checking if ggraph is installed, if not installing it
if(!requireNamespace("ggraph", quietly = TRUE)){
  install.packages("ggraph") 
}
# Checking if ggnewscale is installed, if not installing it
if(!requireNamespace("ggnewscale", quietly = TRUE)){
  install.packages("ggnewscale") 
}

# install IMPRINTS.CETSA from github 
devtools::install_github("nkdailingyun/IMPRINTS.CETSA")
# install IMPRINTS.CETSA.app from github 
devtools::install_github("mgerault/IMPRINTS.CETSA.app")

# if you have any problems during installation, feel free to send an e-mail at:
# marc.antoine.gerault@ki.se

# all data presents in the Data folder are preprocessed normalized data
# if you want the raw data and renalyze them from scratch, follow the link in the paper 
# and the method explained in it

# all scripts present in this repository are enough to reproduce analysis and figures from the paper



