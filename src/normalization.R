# Normalization functions

normalize_expression_data <- function(data) {
  required_components <- c("filtered_data", "metadata")
  missing_components <- setdiff(required_components, names(data))

  if (length(missing_components) > 0L) {
    stop(
      "Missing required data components: ",
      paste(missing_components, collapse = ", ")
    )
  }

  if (!identical(colnames(data$filtered_data), rownames(data$metadata))) {
    stop(
      "The samples in 'filtered_data' and 'metadata' must have the same ",
      "names and order before normalization."
    )
  }

  # Estimate size factors once so that both output matrices are derived from
  # the same DESeq2 object.
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = data$filtered_data,
    colData = data$metadata,
    design = ~1
  )
  dds <- DESeq2::estimateSizeFactors(dds)

  # Keep size-factor-normalized counts for downstream expression plots.
  normalized_counts <- as.data.frame(
    DESeq2::counts(dds, normalized = TRUE),
    check.names = FALSE
  )

  # WGCNA expects samples in rows and genes in columns.
  dds_vst <- DESeq2::vst(dds)
  expression_vst <- t(SummarizedExperiment::assay(dds_vst))

  list(
    expression_vst = expression_vst,
    normalized_counts = normalized_counts
  )
}

# Retain the original helper for compatibility with existing analyses.
vst_normalization <- function(data) {
  normalize_expression_data(data)$expression_vst
}


