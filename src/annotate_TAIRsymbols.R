if (!("org.At.tair.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.At.tair.db", update = FALSE)
}

library(org.At.tair.db)

get_TAIRsymbols <- function(genes) {
  x <- org.At.tairSYMBOL
  mapped_tairs <- mappedkeys(x)
  xx <- as.list(x[mapped_tairs])
  annotation <- vector("character", length(genes))
  for (j in seq_along(genes)) {
    i <- genes[j]
    annotation[j] <- if (!is.null(xx[[i]][1])) xx[[i]][1] else i
  }
  
  return(annotation)
}

TAIR_to_ENTREZ <- function(genes) {
  x <- org.At.tairENTREZID
  mapped_tairs <- mappedkeys(x)
  xx <- as.list(x[mapped_tairs])
  annotation <- vector("character", length(genes))
  for (j in seq_along(genes)) {
    i <- genes[j]
    annotation[j] <- if (!is.null(xx[[i]][1])) xx[[i]][1] else i
  }
  
  return(annotation)
}

ENTREZ_to_TAIR <- function(entrez_ids) {
  library(org.At.tair.db)
  
  # Reverse the mapping from ENTREZID to TAIR
  x <- revmap(org.At.tairENTREZID)
  
  # Get all mapped ENTREZ IDs
  mapped_entrez <- mappedkeys(x)
  
  # Convert to a list for lookup
  xx <- as.list(x[mapped_entrez])
  
  # Initialize a vector to store the results
  annotation <- vector("character", length(entrez_ids))
  
  for (j in seq_along(entrez_ids)) {
    i <- entrez_ids[j]
    annotation[j] <- if (!is.null(xx[[i]][1])) xx[[i]][1] else i
  }
  
  return(annotation)
}