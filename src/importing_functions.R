#functions for data importing

import_counts <- function(filename, gene_length_col, sep) {
  file_path <- file.path("in", filename)
  if (!file.exists(file_path)) {
    stop(paste("Counts file not found:", file_path))
  }
  
  data <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  if (!is.null(gene_length_col)) {
    if (is.numeric(gene_length_col)) {
      if (gene_length_col > ncol(data)) {
        stop("Gene length column index exceeds the number of columns in the data.")
      }
      gene_length <- data[, gene_length_col]
      data <- data[, -gene_length_col]
      sample_names <- colnames(data)
    } 
    else if (is.character(gene_length_col)) {
      if (!(gene_length_col %in% colnames(data))) {
        stop(paste("Gene length column not found:", gene_length_col))
      }
      gene_length <- data[[gene_length_col]]
      data <- data[, !(colnames(data) %in% gene_length_col)]
      sample_names <- colnames(data)
    } else {
      stop("Invalid `gene_length_col`. Must be either numeric or character.")
    }
    
  }
  
  return(list(data = data, 
              gene_length = if (exists("gene_length")) gene_length else NULL,
              samples = sample_names))
}

import_metadata <- function(filename, sep, sample_id) {
  file_path <- file.path("in", filename)
  if (!file.exists(file_path)) {
    stop(paste("Metadata file not found:", file_path))
  }
  
  metadata <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = TRUE)
  
  if (!(sample_id %in% colnames(metadata))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(metadata) <- metadata[[sample_id]]
  return(metadata)
}

import_traits <- function(filename, sep, sample_id) {
  file_path <- file.path("in", filename)
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  if (!(sample_id %in% colnames(data))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(data) <- data[[sample_id]]
  return(data)
}


import_data <- function(sample_names, 
                        counts_filename, gene_length_column = NULL, counts_sep = ",",
                        metadata_filename, metadata_sep = ",", metadata_sample_id = "sample",
                        traits_filename = NULL, traits_sep = ",", traits_sample_id = "sample",){
  
  out <- import_counts(filename = counts_filename, gene_length_col = gene_length_column, sep = counts_sep)
  metadata <- import_metadata(filename = metadata_filename, sep = metadata_sep, sample_id = metadata_sample_id)
  if(metadata$sample != out$samples) {
    stop("metadata sample names do not corespond to data sample names.")
  }
  else{
    out <- c(out, as.list(metadata[,!(metadata_sample_id)]))
  }
  
  if (!(is.null(traits_filename))){
    traits <- import_traits(filename = traits_filename, sep = traits_sep, sample_id = traits_sample_id)
    
    out <- c(out, as.list(traits[, !(traits_sample_id)]))
  }
  
  return(out)
  
}
