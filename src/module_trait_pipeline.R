# Callable implementation of pipeline step 4.

module_trait_network_artifact <- function(session_dir) {
  candidates <- file.path(
    session_dir,
    "out",
    c("3_networks.rds", "3_networks.Rdata", "networks.Rdata")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) NULL else existing[[1L]]
}

read_module_trait_network_result <- function(path) {
  rds_result <- tryCatch(readRDS(path), error = function(error) NULL)
  if (is.list(rds_result) && !is.null(rds_result$network)) {
    return(rds_result)
  }

  network_environment <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    load(path, envir = network_environment),
    error = function(error) character()
  )
  if ("network_result" %in% loaded &&
      is.list(network_environment$network_result) &&
      !is.null(network_environment$network_result$network)) {
    return(network_environment$network_result)
  }
  if ("net" %in% loaded) {
    return(list(
      network = network_environment$net,
      network_kmeans = if ("net_km" %in% loaded) {
        network_environment$net_km
      } else {
        NULL
      }
    ))
  }
  stop("The network artifact has an unsupported format: ", path, call. = FALSE)
}

load_module_trait_inputs <- function(session_dir) {
  prepared_path <- file.path(session_dir, "out", "1_prepared_data.RData")
  network_path <- module_trait_network_artifact(session_dir)
  if (!file.exists(prepared_path)) {
    stop("Prepared data are missing. Complete step 1 first.", call. = FALSE)
  }
  if (is.null(network_path)) {
    stop("Network results are missing. Complete step 3 first.", call. = FALSE)
  }

  prepared_environment <- new.env(parent = emptyenv())
  loaded <- load(prepared_path, envir = prepared_environment)
  if (!"data" %in% loaded || !is.list(prepared_environment$data)) {
    stop("The prepared-data artifact does not contain 'data'.", call. = FALSE)
  }

  list(
    data = prepared_environment$data,
    network_result = read_module_trait_network_result(network_path)
  )
}

observed_trait_levels <- function(values) {
  values <- values[!is.na(values)]
  if (is.factor(values)) {
    return(levels(droplevels(values)))
  }
  unique(as.character(values))
}

module_trait_control_spec <- function(session_dir) {
  inputs <- load_module_trait_inputs(session_dir)
  data <- inputs$data
  metadata <- data$metadata
  traits <- data$traits
  if (!is.data.frame(metadata)) {
    stop("Prepared metadata are missing or invalid.", call. = FALSE)
  }

  categorical <- names(metadata)[vapply(
    metadata,
    function(column) {
      is.factor(column) || is.character(column) || is.logical(column)
    },
    logical(1)
  )]
  categorical_levels <- stats::setNames(
    lapply(metadata[categorical], observed_trait_levels),
    categorical
  )
  categorical_counts <- vapply(categorical_levels, length, integer(1))
  informative_categorical <- categorical[categorical_counts >= 2L]
  skipped_categorical <- categorical[categorical_counts < 2L]

  continuous <- if (is.null(traits)) character() else {
    names(traits)[vapply(traits, is.numeric, logical(1))]
  }
  continuous_counts <- if (length(continuous) == 0L) {
    integer()
  } else {
    vapply(traits[continuous], function(column) {
      length(unique(column[is.finite(column)]))
    }, integer(1))
  }
  informative_continuous <- continuous[continuous_counts >= 2L]
  skipped_continuous <- continuous[continuous_counts < 2L]

  networks <- c(hc = "Hierarchical clustering")
  if (!is.null(inputs$network_result$network_kmeans)) {
    networks <- c(networks, kmeans = "K-means refinement")
  }

  list(
    samples = nrow(metadata),
    metadata = metadata,
    traits = traits,
    categorical = informative_categorical,
    categorical_levels = categorical_levels,
    categorical_counts = categorical_counts,
    skipped_categorical = skipped_categorical,
    continuous = informative_continuous,
    continuous_counts = continuous_counts,
    skipped_continuous = skipped_continuous,
    networks = networks
  )
}

validate_module_trait_network <- function(network, samples, label) {
  if (is.null(network) || !is.matrix(network$MEs) && !is.data.frame(network$MEs)) {
    stop(label, " does not contain a valid eigengene matrix.", call. = FALSE)
  }
  eigengenes <- as.data.frame(network$MEs, check.names = FALSE)
  if (is.null(rownames(eigengenes))) {
    stop(label, " eigengenes do not have sample names.", call. = FALSE)
  }
  if (!setequal(rownames(eigengenes), samples)) {
    stop(label, " eigengenes and prepared data contain different samples.", call. = FALSE)
  }
  eigengenes[samples, , drop = FALSE]
}

module_trait_correlations <- function(eigengenes,
                                      traits,
                                      method = c("pearson", "spearman")) {
  method <- match.arg(method)
  eigengenes <- as.data.frame(eigengenes, check.names = FALSE)
  traits <- as.data.frame(traits, check.names = FALSE)
  correlations <- matrix(
    NA_real_,
    nrow = ncol(eigengenes),
    ncol = ncol(traits),
    dimnames = list(colnames(eigengenes), colnames(traits))
  )
  p_values <- correlations
  sample_counts <- correlations

  for (module_index in seq_len(ncol(eigengenes))) {
    module_values <- eigengenes[[module_index]]
    for (trait_index in seq_len(ncol(traits))) {
      trait_values <- traits[[trait_index]]
      complete <- stats::complete.cases(module_values, trait_values) &
        is.finite(module_values) & is.finite(trait_values)
      sample_counts[module_index, trait_index] <- sum(complete)
      module_sd <- stats::sd(module_values[complete])
      trait_sd <- stats::sd(trait_values[complete])
      if (sum(complete) < 3L || !is.finite(module_sd) ||
          !is.finite(trait_sd) || module_sd == 0 || trait_sd == 0) {
        next
      }
      test <- suppressWarnings(stats::cor.test(
        module_values[complete],
        trait_values[complete],
        method = method,
        exact = FALSE
      ))
      correlations[module_index, trait_index] <- unname(test$estimate)
      p_values[module_index, trait_index] <- test$p.value
    }
  }

  list(correlations = correlations, p_values = p_values, n = sample_counts)
}

module_trait_heatmap <- function(result, title) {
  correlations <- result$correlations
  p_values <- result$p_values
  plot_data <- data.frame(
    module = rep(rownames(correlations), times = ncol(correlations)),
    trait = rep(colnames(correlations), each = nrow(correlations)),
    correlation = as.vector(correlations),
    p_value = as.vector(p_values),
    stringsAsFactors = FALSE
  )
  plot_data$label <- ifelse(
    is.na(plot_data$correlation),
    "NA",
    paste0(
      sprintf("%.2f", plot_data$correlation),
      "\n",
      ifelse(is.na(plot_data$p_value), "p=NA", paste0("p=", format.pval(
        plot_data$p_value, digits = 2, eps = 0.001
      )))
    )
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = trait, y = module, fill = correlation)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6, lineheight = 0.9) +
    ggplot2::scale_fill_gradient2(
      low = "#2563eb",
      mid = "white",
      high = "#dc2626",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "#e5e7eb",
      name = "Correlation"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

write_module_trait_matrix <- function(matrix, path) {
  output <- data.frame(module = rownames(matrix), matrix, check.names = FALSE)
  utils::write.table(
    output,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  normalizePath(path, mustWork = TRUE)
}

run_module_trait_analysis <- function(session_dir,
                                      continuous_traits = character(),
                                      categorical_variables = character(),
                                      reference_levels = NULL,
                                      include_pairwise = TRUE,
                                      include_vs_all = TRUE,
                                      min_count = 1,
                                      drop_redundant = TRUE,
                                      selected_contrasts = NULL,
                                      network_types = c("hc", "kmeans"),
                                      continuous_method = "pearson",
                                      create_plots = TRUE,
                                      save_artifacts = TRUE,
                                      prepared_data = NULL,
                                      network_result = NULL,
                                      verbose = TRUE) {
  if (is.null(prepared_data) || is.null(network_result)) {
    inputs <- load_module_trait_inputs(session_dir)
    prepared_data <- prepared_data %||moduletrait% inputs$data
    network_result <- network_result %||moduletrait% inputs$network_result
  }
  continuous_traits <- unique(continuous_traits)
  categorical_variables <- unique(categorical_variables)
  network_types <- unique(network_types)
  if (length(continuous_traits) == 0L && length(categorical_variables) == 0L) {
    stop("Select at least one continuous or categorical trait.", call. = FALSE)
  }
  if (length(network_types) == 0L ||
      !all(network_types %in% c("hc", "kmeans"))) {
    stop("Select at least one available network.", call. = FALSE)
  }
  if (!continuous_method %in% c("pearson", "spearman")) {
    stop("Continuous correlation must be 'pearson' or 'spearman'.", call. = FALSE)
  }
  if (!is.numeric(min_count) || length(min_count) != 1L ||
      is.na(min_count) || min_count < 1 || min_count %% 1 != 0) {
    stop("Minimum observations per contrast must be a positive integer.", call. = FALSE)
  }

  metadata <- prepared_data$metadata
  traits <- prepared_data$traits
  samples <- rownames(metadata)
  if (is.null(samples) || anyDuplicated(samples)) {
    stop("Prepared metadata must have unique sample row names.", call. = FALSE)
  }
  if (length(continuous_traits) > 0L) {
    if (is.null(traits) || !all(continuous_traits %in% names(traits))) {
      stop("One or more selected continuous traits are unavailable.", call. = FALSE)
    }
    if (!all(vapply(traits[continuous_traits], is.numeric, logical(1)))) {
      stop("Selected continuous traits must be numeric.", call. = FALSE)
    }
    if (!setequal(rownames(traits), samples)) {
      stop("Prepared traits and metadata contain different samples.", call. = FALSE)
    }
    traits <- traits[samples, continuous_traits, drop = FALSE]
  }

  binary_matrix <- NULL
  skipped_categorical <- character()
  dropped_redundant <- character()
  if (length(categorical_variables) > 0L) {
    binary_matrix <- binary_contrast_matrix(
      data = metadata,
      variables = categorical_variables,
      reference_levels = reference_levels,
      vsAll = isTRUE(include_vs_all),
      includePairwise = isTRUE(include_pairwise),
      minCount = as.integer(min_count),
      drop_redundant = isTRUE(drop_redundant),
      verbose = FALSE
    )
    skipped_categorical <- attr(binary_matrix, "skipped_single_level_variables")
    dropped_redundant <- attr(binary_matrix, "dropped_redundant_contrasts")
    binary_matrix <- select_binary_contrasts(
      binary_matrix,
      include = selected_contrasts,
      verbose = FALSE
    )
  }

  networks <- list(hc = network_result$network)
  if (!is.null(network_result$network_kmeans)) {
    networks$kmeans <- network_result$network_kmeans
  }
  unavailable_networks <- setdiff(network_types, names(networks))
  if (length(unavailable_networks) > 0L) {
    stop(
      "Selected networks are unavailable: ",
      paste(unavailable_networks, collapse = ", "),
      call. = FALSE
    )
  }
  network_labels <- c(hc = "Hierarchical clustering", kmeans = "K-means refinement")

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(binary_contrasts = NULL, report = NULL, tables = list(), plots = list())

  if (save_artifacts && !is.null(binary_matrix)) {
    binary_output <- data.frame(
      sample = rownames(binary_matrix),
      binary_matrix,
      check.names = FALSE
    )
    files$binary_contrasts <- file.path(
      output_dir, "module_trait_binary_contrasts.tsv"
    )
    utils::write.table(
      binary_output,
      files$binary_contrasts,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  results <- list()
  for (network_name in network_types) {
    eigengenes <- validate_module_trait_network(
      networks[[network_name]], samples, network_labels[[network_name]]
    )
    network_results <- list(continuous = NULL, categorical = NULL)

    if (length(continuous_traits) > 0L) {
      continuous_result <- module_trait_correlations(
        eigengenes,
        traits,
        method = continuous_method
      )
      continuous_result$plot <- module_trait_heatmap(
        continuous_result,
        paste(network_labels[[network_name]], "- continuous traits")
      )
      network_results$continuous <- continuous_result
    }
    if (!is.null(binary_matrix)) {
      categorical_result <- module_trait_correlations(
        eigengenes,
        binary_matrix,
        method = "pearson"
      )
      categorical_result$plot <- module_trait_heatmap(
        categorical_result,
        paste(network_labels[[network_name]], "- categorical contrasts")
      )
      network_results$categorical <- categorical_result
    }
    results[[network_name]] <- network_results

    if (save_artifacts) {
      for (analysis_name in c("continuous", "categorical")) {
        current <- network_results[[analysis_name]]
        if (is.null(current)) next
        prefix <- paste("4", network_name, analysis_name, sep = "_")
        correlation_path <- file.path(
          output_dir, paste0(prefix, "_correlations.tsv")
        )
        p_value_path <- file.path(output_dir, paste0(prefix, "_p_values.tsv"))
        files$tables[[paste0(network_name, "_", analysis_name)]] <- list(
          correlations = write_module_trait_matrix(
            current$correlations, correlation_path
          ),
          p_values = write_module_trait_matrix(current$p_values, p_value_path)
        )
        if (create_plots) {
          plot_path <- if (identical(analysis_name, "categorical")) {
            file.path(
              plots_dir,
              if (identical(network_name, "hc")) {
                "EMhc-condition_corr.pdf"
              } else {
                "EMkm-condition_corr.pdf"
              }
            )
          } else {
            file.path(
              plots_dir,
              paste0("4_", network_name, "_continuous_traits_correlation.pdf")
            )
          }
          ggplot2::ggsave(
            plot_path,
            current$plot,
            width = max(8, min(18, ncol(current$correlations) * 1.25)),
            height = max(6, min(24, nrow(current$correlations) * 0.38 + 2)),
            bg = "white",
            limitsize = FALSE
          )
          files$plots[[paste0(network_name, "_", analysis_name)]] <- plot_path
        }
      }
    }
  }

  parameters <- list(
    continuous_traits = continuous_traits,
    categorical_variables = categorical_variables,
    reference_levels = reference_levels,
    include_pairwise = isTRUE(include_pairwise),
    include_vs_all = isTRUE(include_vs_all),
    min_count = as.integer(min_count),
    drop_redundant = isTRUE(drop_redundant),
    selected_contrasts = if (is.null(binary_matrix)) character() else {
      colnames(binary_matrix)
    },
    network_types = network_types,
    continuous_method = continuous_method
  )
  report <- list(
    parameters = parameters,
    samples = length(samples),
    modules = stats::setNames(
      vapply(results, function(result) {
        current <- result$continuous %||moduletrait% result$categorical
        nrow(current$correlations)
      }, integer(1)),
      network_types
    ),
    contrasts = if (is.null(binary_matrix)) 0L else ncol(binary_matrix),
    skipped_categorical = skipped_categorical,
    dropped_redundant = dropped_redundant,
    results = results,
    files = files,
    completed_at = Sys.time()
  )
  class(report) <- c("wgcnamp_module_trait_report", class(report))
  if (save_artifacts) {
    files$report <- file.path(output_dir, "4_module_trait_report.rds")
    report$files$report <- files$report
    saveRDS(report, files$report)
  }
  if (verbose) {
    message(
      "Module–trait analysis complete for ", length(network_types),
      " network(s), ", length(continuous_traits), " continuous trait(s), and ",
      report$contrasts, " categorical contrast(s)."
    )
  }
  report
}

load_module_trait_report <- function(session_dir) {
  path <- file.path(session_dir, "out", "4_module_trait_report.rds")
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(error) NULL)
}

print.wgcnamp_module_trait_report <- function(x, ...) {
  cat("WGCNAmp module–trait analysis\n")
  cat("  Samples:    ", x$samples, "\n", sep = "")
  cat("  Networks:   ", length(x$modules), "\n", sep = "")
  cat("  Contrasts:  ", x$contrasts, "\n", sep = "")
  invisible(x)
}

`%||moduletrait%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}
