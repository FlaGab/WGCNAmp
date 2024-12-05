#functions for data importing 

import_counts <- function(file_path, gene_length_col = NULL, sep = ",") {
  if (!file.exists(file_path)) {
    stop(paste("Counts file not found:", file_path))
  }
  
  data <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  gene_length <- NULL
  if (!is.null(gene_length_col)) {
    if (is.numeric(gene_length_col)) {
      if (gene_length_col > ncol(data)) {
        stop("Gene length column index exceeds the number of columns in the data.")
      }
      gene_length <- data[, gene_length_col]
      data <- data[, -gene_length_col, drop = FALSE]
    } else if (is.character(gene_length_col)) {
      if (!(gene_length_col %in% colnames(data))) {
        stop(paste("Gene length column not found:", gene_length_col))
      }
      gene_length <- data[[gene_length_col]]
      data <- data[, !(colnames(data) %in% gene_length_col), drop = FALSE]
    } else {
      stop("Invalid `gene_length_col`. Must be either numeric or character.")
    }
  }
  
  return(list(data = data, gene_length = gene_length, samples = colnames(data)))
}

import_metadata <- function(file_path, sep = ",", sample_id = "sample") {
  if (!file.exists(file_path)) {
    stop(paste("Metadata file not found:", file_path))
  }
  
  metadata <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = TRUE)
  
  if (!(sample_id %in% colnames(metadata))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(metadata) <- metadata[[sample_id]]
  metadata <- metadata[, colnames(metadata)!=sample_id]
  return(metadata)
}

import_traits <- function(file_path, sep = ",", sample_id = "sample") {
  if (!file.exists(file_path)) {
    stop(paste("Traits file not found:", file_path))
  }
  
  traits <- read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  if (!(sample_id %in% colnames(traits))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(traits) <- traits[[sample_id]]
  traits <- traits %>% select(-all_of(sample_id))
  return(traits)
}

import_data <- function(counts_file, gene_length_column = NULL, counts_sep = ",", 
                        metadata_file, metadata_sep = ",", metadata_sample_id = "sample", 
                        traits_file = NULL, traits_sep = ",", traits_sample_id = "sample") {
  # Import counts
  counts <- import_counts(file_path = counts_file, gene_length_col = gene_length_column, sep = counts_sep)
  
  # Import metadata
  metadata <- import_metadata(file_path = metadata_file, sep = metadata_sep, sample_id = metadata_sample_id)
  
  # Check sample correspondence
  if (!all(colnames(counts$data) %in% rownames(metadata))) {
    stop("Mismatch between metadata sample names and counts sample names.")
  }
  
  # Optional: Import traits
  traits <- NULL
  if (!is.null(traits_file)) {
    traits <- import_traits(file_path = traits_file, sep = traits_sep, sample_id = traits_sample_id)
    if (!all(rownames(metadata) %in% rownames(traits))) {
      stop("Mismatch between metadata sample names and traits sample names.")
    }
  }
  
  # Combine and return
  return(list(
    counts = counts$data,
    gene_length = counts$gene_length,
    metadata = metadata[,!(colnames(metadata)==metadata_sample_id)],
    traits = if (!is.null(traits)) traits else NULL
  ))
}