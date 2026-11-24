#functions for data importing 

import_counts <- function(file_path, gene_length_col = NULL, sep = ",") {
  if (!file.exists(file_path)) {
    stop(paste("Counts file not found:", file_path))
  }
  
  data <- read.table(
    file_path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Standard CSV writers commonly store gene identifiers in the first column
  # (often with an empty header) instead of encoding them as row names. Accept
  # that layout automatically while preserving the legacy row-name layout.
  default_row_names <- identical(rownames(data), as.character(seq_len(nrow(data))))
  if (default_row_names && ncol(data) > 0L && !is.numeric(data[[1L]])) {
    gene_ids <- as.character(data[[1L]])
    if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
      stop("The gene-identifier column contains missing or empty values.")
    }
    if (anyDuplicated(gene_ids)) {
      stop("Counts data contains duplicated gene identifiers.")
    }
    rownames(data) <- gene_ids
    data <- data[, -1L, drop = FALSE]
  }
  
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

  if (anyDuplicated(rownames(data))) {
    stop("Counts data contains duplicated gene identifiers.")
  }
  if (anyDuplicated(colnames(data))) {
    stop("Counts data contains duplicated sample identifiers.")
  }
  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("Every expression column in the counts file must be numeric.")
  }

  if (!is.null(gene_length)) {
    names(gene_length) <- rownames(data)
  }
  
  return(list(data = data, gene_length = gene_length, samples = colnames(data)))
}

import_metadata <- function(file_path, sep = ",", sample_id = "sample") {
  if (!file.exists(file_path)) {
    stop(paste("Metadata file not found:", file_path))
  }
  
  metadata <- read.table(
    file_path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = TRUE,
    check.names = FALSE
  )
  
  if (!(sample_id %in% colnames(metadata))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(metadata) <- metadata[[sample_id]]
  metadata <- metadata[, colnames(metadata) != sample_id, drop = FALSE]
  if (anyDuplicated(rownames(metadata))) {
    stop("Metadata contains duplicated sample identifiers.")
  }
  return(metadata)
}

import_traits <- function(file_path, sep = ",", sample_id = "sample") {
  if (!file.exists(file_path)) {
    stop(paste("Traits file not found:", file_path))
  }
  
  traits <- read.table(
    file_path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!(sample_id %in% colnames(traits))) {
    stop(paste("Sample ID column not found:", sample_id))
  }
  
  rownames(traits) <- traits[[sample_id]]
  traits <- traits[, setdiff(colnames(traits), sample_id), drop = FALSE]
  if (anyDuplicated(rownames(traits))) {
    stop("Traits data contains duplicated sample identifiers.")
  }
  return(traits)
}

import_data <- function(counts_file, gene_length_column = NULL, counts_sep = ",", 
                        metadata_file, metadata_sep = ",", metadata_sample_id = "sample", 
                        traits_file = NULL, traits_sep = ",", traits_sample_id = "sample") {
  # Import counts
  counts <- import_counts(file_path = counts_file, gene_length_col = gene_length_column, sep = counts_sep)
  
  # Import metadata
  metadata <- import_metadata(file_path = metadata_file, sep = metadata_sep, sample_id = metadata_sample_id)
  
  # Optional: Import traits before validating and aligning all sample tables.
  traits <- NULL
  if (!is.null(traits_file)) {
    traits <- import_traits(
      file_path = traits_file,
      sep = traits_sep,
      sample_id = traits_sample_id
    )
  }

  # Check sample correspondence and enforce count-matrix sample order.
  sample_ids <- colnames(counts$data)

  if (!setequal(sample_ids, rownames(metadata))) {
    stop("Counts and metadata contain different sample identifiers.")
  }

  metadata <- metadata[sample_ids, , drop = FALSE]

  if (!is.null(traits)) {
    if (!setequal(sample_ids, rownames(traits))) {
      stop("Counts and traits contain different sample identifiers.")
    }
    traits <- traits[sample_ids, , drop = FALSE]
  }

  stopifnot(
    identical(sample_ids, rownames(metadata)),
    is.null(traits) || identical(sample_ids, rownames(traits))
  )
  # Combine and return
  return(list(
    counts = counts$data,
    gene_length = counts$gene_length,
    metadata = metadata,
    traits = traits
  ))
}
