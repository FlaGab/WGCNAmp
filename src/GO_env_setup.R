#GO enrichmant with clusterprofiler evironment setup

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("org.At.tair.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.At.tair.db", update = FALSE)
}

if (!("GO.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GO.db", update = FALSE)
}

if (!("openxlsx" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("openxlsx")
}

if (!("enrichplot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("enrichplot", update = FALSE)
}

library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(GO.db)
library(openxlsx)
library(enrichplot)
