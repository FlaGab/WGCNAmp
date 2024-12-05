#Environment setup for WGCNA
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("devtools" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("devtools")
}

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ComplexHeatmap")
}

if (!("dplyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("dplyr")
}

if (!("tidyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("tidyr")
}

if (!("readr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("readr")
}

if (!("limma" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("limma")
}

if (!("sva" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("sva")
}

if (!("CorLevelPlot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github("kevinblighe/CorLevelPlot")
}
  
if (!requireNamespace(c("graph", "RBGL", "topGO", "GOSim"), quietly = TRUE)){
  BiocManager::install(c("graph", "RBGL", "topGO", "GOSim"))
}

if (!("CoExpNets" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github('juanbot/CoExpNets', force = T)
}

if (!("ggdendro" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggdendro")
}

if (!("ggrepel" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggrepel")
}

if (!("dendextend" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("dendextend")
}

library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
library(sva)
library(gridExtra)
library(CorLevelPlot)
library(CoExpNets)
library(ggdendro)
library(ggrepel)
library(dendextend)
source("src/make_module_heatmap.R")
source("src/export_to_cytoscape.R")
source("src/plot_genes_in_modules.R")
source("src/col2num.R")
source("src/annotate_TAIRsymbols.R")

source("src/importing_functions_v2.R")
source("src/quality_check.R")
source("src/filtering_functions.R")
source("src/normalization.R")
