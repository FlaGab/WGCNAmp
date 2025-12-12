#environment setup for running module 7_visualization_genes
library(DESeq2)
library(ggplot2)
library(pheatmap)

load(paste0(session_name, "/out/1_prepared_data.Rdata"))
load(paste0(session_name, "/out/3_networks.Rdata"))
source(paste0(session_name, "/in/", session_name,"_parameters.R"))