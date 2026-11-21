# Callable implementation of pipeline step 1.

save_pdf_output <- function(path, draw, width = 8, height = 6) {
  grDevices::pdf(
    file = path,
    width = width,
    height = height,
    bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

run_data_preparation <- function(session_dir,
                                 counts_file,
                                 metadata_file,
                                 gene_length_column = NULL,
                                 traits_file = NULL,
                                 excluded_samples = character(),
                                 filtering_method = "2",
                                 n = 10,
                                 m = 3,
                                 p = 0.75,
                                 counts_sep = ",",
                                 metadata_sep = ",",
                                 metadata_sample_id = "sample",
                                 traits_sep = ",",
                                 traits_sample_id = "sample",
                                 qc_color_column = NULL,
                                 metadata_transform = NULL,
                                 create_plots = TRUE,
                                 save_artifact = TRUE,
                                 normalization_function = normalize_expression_data,
                                 verbose = TRUE) {
  required_functions <- c(
    "import_data",
    "missing_QC",
    "HC_QC",
    "PCA_QC",
    "exclude_samples",
    "low_reads_filtering"
  )
  helper_environment <- environment()
  missing_functions <- required_functions[!vapply(
    required_functions,
    function(name) {
      exists(
        name,
        envir = helper_environment,
        mode = "function",
        inherits = TRUE
      )
    },
    logical(1)
  )]
  if (length(missing_functions) > 0L) {
    stop(
      "Data-preparation helpers have not been loaded: ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(session_dir) || length(session_dir) != 1L || !nzchar(session_dir)) {
    stop("'session_dir' must be one non-empty path.", call. = FALSE)
  }
  if (!is.character(excluded_samples)) {
    stop("'excluded_samples' must be a character vector.", call. = FALSE)
  }
  if (!is.null(metadata_transform) && !is.function(metadata_transform)) {
    stop("'metadata_transform' must be NULL or a function.", call. = FALSE)
  }
  if (!is.function(normalization_function)) {
    stop("'normalization_function' must be a function.", call. = FALSE)
  }

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message("Importing expression data and sample annotations...")
  }
  data <- import_data(
    counts_file = counts_file,
    gene_length_column = gene_length_column,
    counts_sep = counts_sep,
    metadata_file = metadata_file,
    metadata_sep = metadata_sep,
    metadata_sample_id = metadata_sample_id,
    traits_file = traits_file,
    traits_sep = traits_sep,
    traits_sample_id = traits_sample_id
  )

  input_gene_count <- nrow(data$counts)
  input_sample_ids <- colnames(data$counts)
  input_sample_count <- length(input_sample_ids)

  unknown_exclusions <- setdiff(excluded_samples, input_sample_ids)
  if (length(unknown_exclusions) > 0L) {
    stop(
      "Excluded samples are not present in the input data: ",
      paste(unknown_exclusions, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(metadata_transform)) {
    transformed_metadata <- metadata_transform(data$metadata)
    if (!is.data.frame(transformed_metadata)) {
      stop("'metadata_transform' must return a data frame.", call. = FALSE)
    }
    if (!setequal(rownames(transformed_metadata), input_sample_ids)) {
      stop(
        "The metadata transformation changed the set of sample identifiers.",
        call. = FALSE
      )
    }
    data$metadata <- transformed_metadata[input_sample_ids, , drop = FALSE]
  }

  if (!is.null(qc_color_column) && !qc_color_column %in% colnames(data$metadata)) {
    stop(
      "QC color column '", qc_color_column, "' is not present in metadata.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Checking genes and samples for missing or invalid values...")
  }
  gsg <- missing_QC(data$counts)
  automatically_excluded_samples <- input_sample_ids[!gsg$goodSamples]
  data$data_no_missing <- data$counts[
    gsg$goodGenes,
    gsg$goodSamples,
    drop = FALSE
  ]
  retained_qc_samples <- colnames(data$data_no_missing)
  data$metadata <- data$metadata[retained_qc_samples, , drop = FALSE]
  if (!is.null(data$traits)) {
    data$traits <- data$traits[retained_qc_samples, , drop = FALSE]
  }
  if (!is.null(data$gene_length)) {
    data$gene_length <- data$gene_length[rownames(data$data_no_missing)]
  }

  if (ncol(data$data_no_missing) < 2L) {
    stop("Fewer than two samples remain after missing-data QC.", call. = FALSE)
  }

  clustering_path <- file.path(plots_dir, "clustering_QC.pdf")
  pca_path <- file.path(plots_dir, "PCA_QC.pdf")
  pca_plot <- PCA_QC(
    data$data_no_missing,
    metadata = data$metadata,
    color_col = qc_color_column
  )

  if (create_plots) {
    hc_tree <- save_pdf_output(clustering_path, function() {
      HC_QC(
        data$data_no_missing,
        metadata = data$metadata,
        color_col = qc_color_column
      )
    })
    save_pdf_output(pca_path, function() {
      print(pca_plot)
      invisible(pca_plot)
    })
  } else {
    hc_tree <- stats::hclust(stats::dist(t(data$data_no_missing)), method = "average")
    clustering_path <- NULL
    pca_path <- NULL
  }

  manually_excluded_samples <- setdiff(
    unique(excluded_samples),
    automatically_excluded_samples
  )
  if (verbose && length(manually_excluded_samples) > 0L) {
    message(
      "Excluding configured samples: ",
      paste(manually_excluded_samples, collapse = ", ")
    )
  }
  data <- exclude_samples(data, samples = manually_excluded_samples)
  if (ncol(data$counts) < 2L) {
    stop("Fewer than two samples remain after configured exclusions.", call. = FALSE)
  }

  if (verbose) {
    message("Filtering low-expression genes...")
  }
  data$filtered_data <- low_reads_filtering(
    data$counts,
    method = filtering_method,
    n = n,
    m = m,
    p = p
  )
  if (nrow(data$filtered_data) == 0L) {
    stop("Expression filtering removed every gene.", call. = FALSE)
  }
  if (!is.null(data$gene_length)) {
    data$gene_length <- data$gene_length[rownames(data$filtered_data)]
  }

  if (verbose) {
    message("Normalizing expression data...")
  }
  normalization_result <- normalization_function(data)
  required_normalization_outputs <- c("expression_vst", "normalized_counts")
  if (!is.list(normalization_result) ||
      !all(required_normalization_outputs %in% names(normalization_result))) {
    stop(
      "The normalization function must return 'expression_vst' and ",
      "'normalized_counts'.",
      call. = FALSE
    )
  }

  expected_genes <- rownames(data$filtered_data)
  expected_samples <- rownames(data$metadata)
  expression_vst <- normalization_result$expression_vst
  normalized_counts <- normalization_result$normalized_counts

  if (!is.matrix(expression_vst) && !is.data.frame(expression_vst)) {
    stop("'expression_vst' must be a matrix or data frame.", call. = FALSE)
  }
  if (!is.matrix(normalized_counts) && !is.data.frame(normalized_counts)) {
    stop("'normalized_counts' must be a matrix or data frame.", call. = FALSE)
  }
  if (!setequal(rownames(expression_vst), expected_samples) ||
      !setequal(colnames(expression_vst), expected_genes)) {
    stop(
      "'expression_vst' must contain samples in rows and retained genes in columns.",
      call. = FALSE
    )
  }
  if (!setequal(rownames(normalized_counts), expected_genes) ||
      !setequal(colnames(normalized_counts), expected_samples)) {
    stop(
      "'normalized_counts' must contain retained genes in rows and samples in columns.",
      call. = FALSE
    )
  }

  data$exprData <- as.matrix(
    expression_vst[expected_samples, expected_genes, drop = FALSE]
  )
  data$normalized_counts <- as.data.frame(
    normalized_counts[expected_genes, expected_samples, drop = FALSE],
    check.names = FALSE
  )

  filtered_gene_count <- nrow(data$filtered_data)
  retained_sample_ids <- colnames(data$filtered_data)
  data$counts <- NULL
  data$filtered_data <- NULL

  artifact_path <- file.path(output_dir, "1_prepared_data.RData")
  if (save_artifact) {
    save(data, file = artifact_path)
  } else {
    artifact_path <- NULL
  }

  summary <- list(
    input_genes = input_gene_count,
    genes_after_missing_qc = sum(gsg$goodGenes),
    retained_genes = filtered_gene_count,
    removed_genes = input_gene_count - filtered_gene_count,
    input_samples = input_sample_count,
    retained_samples = length(retained_sample_ids),
    automatic_sample_exclusions = automatically_excluded_samples,
    requested_sample_exclusions = unique(excluded_samples),
    configured_sample_exclusions = manually_excluded_samples,
    retained_sample_ids = retained_sample_ids,
    filtering = list(method = filtering_method, n = n, m = m, p = p)
  )

  result <- list(
    data = data,
    qc = list(
      good_samples_genes = gsg,
      clustering = hc_tree,
      pca = pca_plot
    ),
    summary = summary,
    files = list(
      prepared_data = artifact_path,
      clustering_qc = clustering_path,
      pca_qc = pca_path,
      preparation_report = if (save_artifact) {
        file.path(output_dir, "1_preparation_report.rds")
      } else {
        NULL
      }
    ),
    completed_at = Sys.time()
  )
  class(result) <- c("wgcnamp_data_preparation", class(result))

  if (save_artifact) {
    report <- result
    report$data <- NULL
    saveRDS(report, result$files$preparation_report)
  }

  if (verbose) {
    message(
      "Data preparation complete: ",
      filtered_gene_count, " genes and ", length(retained_sample_ids),
      " samples retained."
    )
  }
  result
}

print.wgcnamp_data_preparation <- function(x, ...) {
  cat("WGCNAmp data preparation\n")
  cat("  Input genes:     ", x$summary$input_genes, "\n", sep = "")
  cat("  Retained genes:  ", x$summary$retained_genes, "\n", sep = "")
  cat("  Input samples:   ", x$summary$input_samples, "\n", sep = "")
  cat("  Retained samples:", x$summary$retained_samples, "\n", sep = "")
  if (!is.null(x$files$prepared_data)) {
    cat("  Artifact:        ", x$files$prepared_data, "\n", sep = "")
  }
  invisible(x)
}
