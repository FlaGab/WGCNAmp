# Callable implementation of pipeline step 2.

default_soft_threshold_powers <- function() {
  c(seq.int(1L, 20L), seq.int(22L, 50L, by = 2L))
}

validate_soft_threshold_powers <- function(powers) {
  if (!is.numeric(powers) || length(powers) == 0L || anyNA(powers) ||
      any(!is.finite(powers)) || any(powers <= 0) || any(powers %% 1 != 0)) {
    stop("'powers' must contain one or more positive integers.", call. = FALSE)
  }
  powers <- sort(unique(as.integer(powers)))
  if (length(powers) < 2L) {
    stop("Evaluate at least two candidate powers.", call. = FALSE)
  }
  powers
}

parse_soft_threshold_powers <- function(specification) {
  if (!is.character(specification) || length(specification) != 1L ||
      !nzchar(trimws(specification))) {
    stop("Enter candidate powers, for example '1:20, 22:50:2'.", call. = FALSE)
  }

  tokens <- trimws(strsplit(specification, ",", fixed = TRUE)[[1L]])
  if (any(!nzchar(tokens))) {
    stop("Candidate powers contain an empty entry.", call. = FALSE)
  }

  parse_token <- function(token) {
    if (grepl("^[0-9]+$", token)) {
      return(as.integer(token))
    }
    pieces <- strsplit(token, ":", fixed = TRUE)[[1L]]
    if (!length(pieces) %in% c(2L, 3L) ||
        any(!grepl("^[0-9]+$", pieces))) {
      stop(
        "Invalid candidate-power entry '", token,
        "'. Use values or ranges such as 1:20 and 22:50:2.",
        call. = FALSE
      )
    }
    values <- as.integer(pieces)
    from <- values[[1L]]
    to <- values[[2L]]
    by <- if (length(values) == 3L) values[[3L]] else 1L
    if (from > to || by <= 0L) {
      stop("Invalid ascending range '", token, "'.", call. = FALSE)
    }
    seq.int(from, to, by = by)
  }

  validate_soft_threshold_powers(unlist(lapply(tokens, parse_token), use.names = FALSE))
}

load_prepared_expression <- function(session_dir) {
  artifact <- file.path(session_dir, "out", "1_prepared_data.RData")
  if (!file.exists(artifact)) {
    stop(
      "Prepared data were not found. Complete step 1 before estimating a soft threshold.",
      call. = FALSE
    )
  }

  environment <- new.env(parent = emptyenv())
  loaded_names <- load(artifact, envir = environment)
  if (!"data" %in% loaded_names || is.null(environment$data$exprData)) {
    stop(
      "The step-1 artifact does not contain 'data$exprData'. Re-run data preparation.",
      call. = FALSE
    )
  }
  environment$data$exprData
}

validate_wgcna_expression <- function(expression_data) {
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("Expression data must be a matrix or data frame.", call. = FALSE)
  }
  expression_data <- as.data.frame(expression_data, check.names = FALSE)
  if (nrow(expression_data) < 2L || ncol(expression_data) < 2L) {
    stop("Expression data must contain at least two samples and two genes.", call. = FALSE)
  }
  if (!all(vapply(expression_data, is.numeric, logical(1)))) {
    stop("Every expression column must be numeric.", call. = FALSE)
  }
  if (anyNA(expression_data) || any(!is.finite(as.matrix(expression_data)))) {
    stop("Expression data contain missing or non-finite values.", call. = FALSE)
  }
  expression_data
}

soft_threshold_fit_table <- function(estimate) {
  if (!is.list(estimate) || is.null(estimate$fitIndices)) {
    stop("The threshold estimator did not return 'fitIndices'.", call. = FALSE)
  }
  fit <- as.data.frame(estimate$fitIndices, check.names = FALSE)
  required <- c("Power", "SFT.R.sq", "slope", "mean.k.")
  missing <- setdiff(required, colnames(fit))
  if (length(missing) > 0L) {
    stop(
      "Soft-threshold fit indices are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  fit$model_fit <- -sign(fit$slope) * fit$SFT.R.sq
  fit
}

select_soft_threshold_power <- function(fit,
                                        power_estimate = NA_real_,
                                        selection_mode = c(
                                          "first_above_cutoff",
                                          "automatic",
                                          "manual"
                                        ),
                                        manual_power = NULL,
                                        fit_cutoff = 0.8) {
  selection_mode <- match.arg(selection_mode)
  if (!is.numeric(fit_cutoff) || length(fit_cutoff) != 1L ||
      is.na(fit_cutoff) || fit_cutoff <= 0 || fit_cutoff > 1) {
    stop("'fit_cutoff' must be a number greater than 0 and at most 1.", call. = FALSE)
  }

  evaluated_powers <- fit$Power
  selected <- switch(
    selection_mode,
    first_above_cutoff = {
      candidates <- fit$Power[is.finite(fit$model_fit) & fit$model_fit >= fit_cutoff]
      if (length(candidates) == 0L) {
        stop(
          "No candidate power reaches the requested model-fit target. ",
          "Review the plot, lower the target, extend the candidates, or select a manual power.",
          call. = FALSE
        )
      }
      min(candidates)
    },
    automatic = {
      if (length(power_estimate) != 1L || !is.finite(power_estimate)) {
        stop(
          "WGCNA did not return an automatic power estimate. ",
          "Use the fit target or select a manual power.",
          call. = FALSE
        )
      }
      power_estimate
    },
    manual = {
      if (!is.numeric(manual_power) || length(manual_power) != 1L ||
          !is.finite(manual_power) || manual_power <= 0 || manual_power %% 1 != 0) {
        stop("Enter one positive integer as the manual power.", call. = FALSE)
      }
      if (!manual_power %in% evaluated_powers) {
        stop(
          "The manual power must be included among the evaluated candidate powers.",
          call. = FALSE
        )
      }
      as.integer(manual_power)
    }
  )
  as.numeric(selected)
}

soft_threshold_plots <- function(fit, selected_power, fit_cutoff, network_type) {
  finite_fit <- fit$model_fit[is.finite(fit$model_fit)]
  lower_limit <- if (length(finite_fit) > 0L) min(0, finite_fit) else 0

  scale_independence <- ggplot2::ggplot(
    fit,
    ggplot2::aes(x = Power, y = model_fit, label = Power)
  ) +
    ggplot2::geom_hline(
      yintercept = fit_cutoff,
      color = "#c23b33",
      linewidth = 0.65,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = selected_power,
      color = "#176b87",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(size = 2, color = "#176b87") +
    ggplot2::geom_text(nudge_y = 0.045, size = 3, check_overlap = TRUE) +
    ggplot2::coord_cartesian(ylim = c(lower_limit, 1.05)) +
    ggplot2::labs(
      title = "Scale independence",
      subtitle = paste("Selected power:", selected_power),
      x = "Soft-threshold power",
      y = paste0(network_type, " model fit (signed R²)")
    ) +
    ggplot2::theme_minimal(base_size = 12)

  mean_connectivity <- ggplot2::ggplot(
    fit,
    ggplot2::aes(x = Power, y = mean.k., label = Power)
  ) +
    ggplot2::geom_vline(
      xintercept = selected_power,
      color = "#176b87",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(size = 2, color = "#176b87") +
    ggplot2::geom_text(nudge_y = 0.04 * max(fit$mean.k., na.rm = TRUE),
                       size = 3, check_overlap = TRUE) +
    ggplot2::labs(
      title = "Mean connectivity",
      x = "Soft-threshold power",
      y = "Mean connectivity"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  list(
    scale_independence = scale_independence,
    mean_connectivity = mean_connectivity
  )
}

save_soft_threshold_pdf <- function(path, plots) {
  grDevices::pdf(path, width = 8, height = 9, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  gridExtra::grid.arrange(
    plots$scale_independence,
    plots$mean_connectivity,
    nrow = 2
  )
  invisible(path)
}

run_soft_threshold_selection <- function(session_dir,
                                         network_type = "signed",
                                         powers = default_soft_threshold_powers(),
                                         fit_cutoff = 0.8,
                                         selection_mode = c(
                                           "first_above_cutoff",
                                           "automatic",
                                           "manual"
                                         ),
                                         manual_power = NULL,
                                         expression_data = NULL,
                                         save_artifact = TRUE,
                                         create_plot = TRUE,
                                         threshold_function = WGCNA::pickSoftThreshold,
                                         verbose = TRUE) {
  selection_mode <- match.arg(selection_mode)
  network_types <- c("signed", "unsigned", "signed hybrid")
  if (!network_type %in% network_types) {
    stop(
      "'network_type' must be one of: ", paste(network_types, collapse = ", "),
      call. = FALSE
    )
  }
  powers <- validate_soft_threshold_powers(powers)
  if (!is.function(threshold_function)) {
    stop("'threshold_function' must be a function.", call. = FALSE)
  }
  if (is.null(expression_data)) {
    expression_data <- load_prepared_expression(session_dir)
  }
  expression_data <- validate_wgcna_expression(expression_data)

  if (verbose) {
    message(
      "Evaluating ", length(powers), " candidate powers for a ",
      network_type, " network..."
    )
  }
  estimate <- threshold_function(
    expression_data,
    networkType = network_type,
    powerVector = powers,
    RsquaredCut = fit_cutoff,
    verbose = if (verbose) 5 else 0
  )
  fit <- soft_threshold_fit_table(estimate)
  power_estimate <- estimate$powerEstimate
  if (is.null(power_estimate) || length(power_estimate) == 0L) {
    power_estimate <- NA_real_
  }
  power_threshold <- select_soft_threshold_power(
    fit = fit,
    power_estimate = power_estimate,
    selection_mode = selection_mode,
    manual_power = manual_power,
    fit_cutoff = fit_cutoff
  )
  plots <- soft_threshold_plots(
    fit = fit,
    selected_power = power_threshold,
    fit_cutoff = fit_cutoff,
    network_type = network_type
  )

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  artifact_path <- if (save_artifact) file.path(output_dir, "2_pt.Rdata") else NULL
  report_path <- if (save_artifact) {
    file.path(output_dir, "2_soft_threshold_report.rds")
  } else {
    NULL
  }
  fit_path <- if (save_artifact) {
    file.path(output_dir, "2_soft_threshold_fit.tsv")
  } else {
    NULL
  }
  plot_path <- if (create_plot) file.path(plots_dir, "sft.pdf") else NULL

  if (create_plot) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Install the 'gridExtra' package to save the soft-threshold PDF.", call. = FALSE)
    }
    save_soft_threshold_pdf(plot_path, plots)
  }
  if (save_artifact) {
    save(power_threshold, network_type, file = artifact_path)
    utils::write.table(
      fit,
      file = fit_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  result <- list(
    power_threshold = power_threshold,
    network_type = network_type,
    powers = powers,
    fit_cutoff = fit_cutoff,
    selection_mode = selection_mode,
    automatic_power = power_estimate,
    fit_indices = fit,
    plots = plots,
    files = list(
      power_threshold = artifact_path,
      report = report_path,
      fit_indices = fit_path,
      plot = plot_path
    ),
    completed_at = Sys.time()
  )
  class(result) <- c("wgcnamp_soft_threshold", class(result))

  if (save_artifact) {
    saveRDS(result, report_path)
  }
  if (verbose) {
    message(
      "Soft-threshold selection complete: power ", power_threshold,
      " (", selection_mode, ")."
    )
  }
  result
}

print.wgcnamp_soft_threshold <- function(x, ...) {
  cat("WGCNAmp soft-threshold selection\n")
  cat("  Network type:   ", x$network_type, "\n", sep = "")
  cat("  Selected power: ", x$power_threshold, "\n", sep = "")
  cat("  Selection mode: ", x$selection_mode, "\n", sep = "")
  invisible(x)
}
