#environment setup for step1_data preparation
if (!("BiocManager" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("BiocManager")
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
if (!("tidyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("tidyr")
}

if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("igraph" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("igraph")
}

if (!("CoExpNets" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github('juanbot/CoExpNets')
}

if (!("CorLevelPlot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github("kevinblighe/CorLevelPlot")
}

library(tidyr)
library(dplyr)
library(ggdendro)
library(ggrepel)
library(dendextend)
library(DESeq2)
library(WGCNA)
source("src/importing_functions_v2.R")
source("src/quality_check.R")
source("src/filtering_functions.R")
source("src/normalization.R")

create_session_env <- function(ses = session_name)
if (!(dir.exists(ses))) {
  dir.create(ses)
  dir.create(paste0(ses, "/out"))
  dir.create(paste0(ses, "/out/TOMfiles"))
  dir.create(paste0(ses, "/out/cytoscape"))
  dir.create(paste0(ses, "/out/GO"))
  dir.create(paste0(ses, "/plots"))
  dir.create(paste0(ses, "/in"))
  file.copy("set_parameters.R", file.path(ses, "in", basename(paste0(ses,"_parameters.R"))))
} else {
  warning("Specified session already exists")
}
