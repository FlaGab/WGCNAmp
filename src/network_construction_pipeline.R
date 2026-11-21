# Callable implementation of pipeline step 3.

network_task_state_path <- function(session_dir) {
  file.path(session_dir, "out", "3_network_task_state.rds")
}

write_network_task_state <- function(session_dir, status, ...) {
  dir.create(file.path(session_dir, "out"), recursive = TRUE, showWarnings = FALSE)
  state <- c(
    list(status = status, updated_at = Sys.time()),
    list(...)
  )
  saveRDS(state, network_task_state_path(session_dir))
  invisible(state)
}

read_network_task_state <- function(session_dir) {
  path <- network_task_state_path(session_dir)
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(error) NULL)
}

load_network_construction_inputs <- function(session_dir) {
  prepared_path <- file.path(session_dir, "out", "1_prepared_data.RData")
  threshold_path <- file.path(session_dir, "out", "2_pt.Rdata")
  if (!file.exists(prepared_path)) {
    stop("Complete step 1 before constructing the network.", call. = FALSE)
  }
  if (!file.exists(threshold_path)) {
    stop("Complete step 2 before constructing the network.", call. = FALSE)
  }

  prepared <- new.env(parent = emptyenv())
  load(prepared_path, envir = prepared)
  if (!exists("data", envir = prepared, inherits = FALSE) ||
      is.null(prepared$data$exprData)) {
    stop(
      "The step-1 artifact does not contain 'data$exprData'. Re-run data preparation.",
      call. = FALSE
    )
  }

  threshold <- new.env(parent = emptyenv())
  load(threshold_path, envir = threshold)
  if (!exists("power_threshold", envir = threshold, inherits = FALSE) ||
      !exists("network_type", envir = threshold, inherits = FALSE)) {
    stop(
      "The step-2 artifact must contain 'power_threshold' and 'network_type'.",
      call. = FALSE
    )
  }

  list(
    expression_data = prepared$data$exprData,
    power_threshold = threshold$power_threshold,
    network_type = threshold$network_type
  )
}

validate_network_expression <- function(expression_data) {
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("Expression data must be a matrix or data frame.", call. = FALSE)
  }
  expression_data <- as.data.frame(expression_data, check.names = FALSE)
  if (nrow(expression_data) < 4L || ncol(expression_data) < 2L) {
    stop(
      "Network construction requires at least four samples and two genes.",
      call. = FALSE
    )
  }
  if (!all(vapply(expression_data, is.numeric, logical(1)))) {
    stop("Every expression column must be numeric.", call. = FALSE)
  }
  if (anyNA(expression_data) || any(!is.finite(as.matrix(expression_data)))) {
    stop("Expression data contain missing or non-finite values.", call. = FALSE)
  }
  if (is.null(colnames(expression_data)) || any(!nzchar(colnames(expression_data)))) {
    stop("Every expression column must have a gene identifier.", call. = FALSE)
  }
  if (anyDuplicated(colnames(expression_data))) {
    stop("Expression data contain duplicated gene identifiers.", call. = FALSE)
  }
  expression_data
}

resolve_tom_type <- function(network_type, tom_type = "auto") {
  network_types <- c("signed", "unsigned", "signed hybrid")
  if (!network_type %in% network_types) {
    stop(
      "Unsupported network type '", network_type, "'.",
      call. = FALSE
    )
  }
  if (identical(tom_type, "auto")) {
    return(if (identical(network_type, "unsigned")) "unsigned" else "signed")
  }
  if (!tom_type %in% c("signed", "unsigned")) {
    stop("'tom_type' must be 'auto', 'signed', or 'unsigned'.", call. = FALSE)
  }
  tom_type
}

ensure_wgcna_attached <- function() {
  # Several WGCNA internals resolve the string "cor" through the search path.
  # Attaching WGCNA ensures its weighted-compatible cor() is selected instead
  # of stats::cor(), including inside a fresh background worker.
  if (!"package:WGCNA" %in% search()) {
    suppressPackageStartupMessages(
      base::library("WGCNA", character.only = TRUE)
    )
  }
  invisible(TRUE)
}

validate_network_parameters <- function(method,
                                        power_threshold,
                                        max_block_size,
                                        correlation_type,
                                        min_module_size,
                                        deep_split,
                                        merge_cut_height,
                                        random_seed) {
  method <- match.arg(method, c("blockwise", "stepwise"))
  if (!is.numeric(power_threshold) || length(power_threshold) != 1L ||
      !is.finite(power_threshold) || power_threshold <= 0) {
    stop("The saved soft-threshold power is invalid.", call. = FALSE)
  }
  if (!is.numeric(max_block_size) || length(max_block_size) != 1L ||
      !is.finite(max_block_size) || max_block_size < 2 || max_block_size %% 1 != 0) {
    stop("'max_block_size' must be an integer of at least 2.", call. = FALSE)
  }
  if (!correlation_type %in% c("pearson", "bicor")) {
    stop("'correlation_type' must be 'pearson' or 'bicor'.", call. = FALSE)
  }
  if (!is.numeric(min_module_size) || length(min_module_size) != 1L ||
      !is.finite(min_module_size) || min_module_size < 2 ||
      min_module_size %% 1 != 0) {
    stop("'min_module_size' must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(deep_split) || length(deep_split) != 1L ||
      !deep_split %in% 0:4) {
    stop("'deep_split' must be an integer from 0 to 4.", call. = FALSE)
  }
  if (!is.numeric(merge_cut_height) || length(merge_cut_height) != 1L ||
      !is.finite(merge_cut_height) || merge_cut_height <= 0 ||
      merge_cut_height >= 1) {
    stop("'merge_cut_height' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(random_seed) || length(random_seed) != 1L ||
      !is.finite(random_seed) || random_seed %% 1 != 0) {
    stop("'random_seed' must be one integer.", call. = FALSE)
  }
  method
}

build_blockwise_network <- function(expression_data,
                                    session_dir,
                                    power_threshold,
                                    network_type,
                                    tom_type,
                                    max_block_size,
                                    correlation_type,
                                    min_module_size,
                                    deep_split,
                                    merge_cut_height,
                                    random_seed,
                                    verbose) {
  ensure_wgcna_attached()
  tom_dir <- file.path(session_dir, "out", "TOMfiles")
  dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
  tom_base <- file.path(tom_dir, "TOM-bw")
  network <- WGCNA::blockwiseModules(
    expression_data,
    maxBlockSize = as.integer(max_block_size),
    corType = correlation_type,
    minModuleSize = as.integer(min_module_size),
    deepSplit = as.integer(deep_split),
    networkType = network_type,
    TOMType = tom_type,
    power = power_threshold,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = tom_base,
    randomSeed = as.integer(random_seed),
    useCorOptionsThroughout = FALSE,
    verbose = if (verbose) 3 else 0
  )
  tom_files <- list.files(
    tom_dir,
    pattern = "^TOM-bw.*[.]RData$",
    full.names = TRUE
  )
  list(
    network = network,
    tom = list(
      method = "blockwise",
      type = tom_type,
      files = normalizePath(tom_files, mustWork = FALSE),
      format = "WGCNA blockwise RData"
    )
  )
}

build_stepwise_network <- function(expression_data,
                                   session_dir,
                                   power_threshold,
                                   network_type,
                                   tom_type,
                                   correlation_type,
                                   min_module_size,
                                   deep_split,
                                   merge_cut_height,
                                   save_adjacency,
                                   verbose) {
  ensure_wgcna_attached()
  tom_dir <- file.path(session_dir, "out", "TOMfiles")
  dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
  cor_function <- if (identical(correlation_type, "bicor")) {
    WGCNA::bicor
  } else {
    WGCNA::cor
  }
  adjacency <- WGCNA::adjacency(
    expression_data,
    type = network_type,
    power = power_threshold,
    corFnc = cor_function,
    corOptions = list(use = "p")
  )
  TOM <- WGCNA::TOMsimilarity(
    adjacency,
    TOMType = tom_type,
    verbose = if (verbose) 2 else 0
  )

  tom_rds <- file.path(tom_dir, "TOM-sbs.rds")
  tom_legacy <- file.path(tom_dir, "TOM-sbs.Rdata")
  saveRDS(TOM, tom_rds)
  save(TOM, file = tom_legacy)
  adjacency_path <- NULL
  if (isTRUE(save_adjacency)) {
    adjacency_path <- file.path(session_dir, "out", "adjacency_matrix.rds")
    saveRDS(adjacency, adjacency_path)
  }
  rm(adjacency)

  dendrogram <- stats::hclust(stats::as.dist(1 - TOM), method = "average")
  unmerged_labels <- dynamicTreeCut::cutreeDynamic(
    dendro = dendrogram,
    distM = 1 - TOM,
    method = "hybrid",
    deepSplit = as.integer(deep_split),
    pamRespectsDendro = TRUE,
    minClusterSize = as.integer(min_module_size),
    verbose = if (verbose) 2 else 0
  )
  unmerged_colors <- WGCNA::labels2colors(unmerged_labels)
  names(unmerged_colors) <- colnames(expression_data)
  correlation_function <- if (identical(correlation_type, "bicor")) {
    WGCNA::bicor
  } else {
    stats::cor
  }
  proper_modules <- setdiff(unique(unmerged_colors), "grey")
  if (length(proper_modules) < 2L) {
    eigengenes <- WGCNA::moduleEigengenes(
      expression_data,
      colors = unmerged_colors
    )$eigengenes
    network <- list(colors = unmerged_colors, MEs = eigengenes)
  } else {
    network <- WGCNA::mergeCloseModules(
      exprData = expression_data,
      colors = unmerged_colors,
      corFnc = correlation_function,
      corOptions = list(use = "p"),
      cutHeight = merge_cut_height,
      verbose = if (verbose) 2 else 0
    )
  }
  if (is.null(network$colors) && !is.null(network$newColors)) {
    network$colors <- network$newColors
  }
  if (is.null(network$MEs) && !is.null(network$newMEs)) {
    network$MEs <- network$newMEs
  }
  names(network$colors) <- colnames(expression_data)
  network$unmergedColors <- unmerged_colors
  network$dendrograms <- list(dendrogram)
  network$blockGenes <- list(seq_len(ncol(expression_data)))
  network$TOMFiles <- tom_legacy

  list(
    network = network,
    tom = list(
      method = "stepwise",
      type = tom_type,
      files = normalizePath(c(tom_rds, tom_legacy), mustWork = TRUE),
      format = c("RDS matrix", "legacy RData variable TOM"),
      adjacency = if (is.null(adjacency_path)) NULL else {
        normalizePath(adjacency_path, mustWork = TRUE)
      }
    )
  )
}

apply_kmeans_refinement <- function(network, expression_data, random_seed) {
  if (!requireNamespace("CoExpNets", quietly = TRUE)) {
    stop(
      "Install the 'CoExpNets' package or disable k-means refinement.",
      call. = FALSE
    )
  }
  set.seed(as.integer(random_seed))
  input_network <- network
  if ("colors" %in% names(input_network)) {
    names(input_network)[names(input_network) == "colors"] <- "moduleColors"
  }
  refined <- CoExpNets::applyKMeans(
    net.file = input_network,
    expr.data = expression_data
  )
  if ("moduleColors" %in% names(refined)) {
    names(refined)[names(refined) == "moduleColors"] <- "colors"
  }
  refined
}

network_module_counts <- function(network, network_label) {
  colors <- network$colors
  if (is.null(colors)) {
    stop("The constructed network does not contain module colors.", call. = FALSE)
  }
  counts <- sort(table(colors), decreasing = TRUE)
  data.frame(
    network = network_label,
    module = names(counts),
    genes = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

network_module_count_plot <- function(counts, title) {
  ggplot2::ggplot(
    counts,
    ggplot2::aes(
      x = stats::reorder(module, genes),
      y = genes,
      fill = module
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = title, x = "Module", y = "Genes") +
    ggplot2::theme_minimal(base_size = 12)
}

save_network_dendrograms <- function(network, session_dir, label) {
  files <- character()
  if (is.null(network$dendrograms) || is.null(network$blockGenes)) {
    return(files)
  }
  for (index in seq_along(network$dendrograms)) {
    file <- file.path(
      session_dir,
      "plots",
      paste0("3_", label, "_modules_dendrogram_block_", index, ".pdf")
    )
    grDevices::pdf(file, width = 9, height = 6, bg = "white")
    device <- grDevices::dev.cur()
    plotted <- tryCatch(
      {
        block_genes <- network$blockGenes[[index]]
        WGCNA::plotDendroAndColors(
          network$dendrograms[[index]],
          network$colors[block_genes],
          paste(label, "modules"),
          dendroLabels = FALSE,
          hang = 0.03,
          addGuide = TRUE,
          guideHang = 0.05
        )
        TRUE
      },
      error = function(error) {
        message(
          "Could not render dendrogram block ", index, ": ",
          conditionMessage(error)
        )
        FALSE
      },
      finally = {
        open_devices <- grDevices::dev.list()
        if (!is.null(open_devices) && device %in% open_devices) {
          grDevices::dev.off(device)
        }
      }
    )
    if (plotted) {
      files <- c(files, file)
    } else if (file.exists(file)) {
      unlink(file)
    }
  }
  normalizePath(files, mustWork = TRUE)
}

save_eigengene_network_plot <- function(network, path, label) {
  if (is.null(network$MEs) || ncol(network$MEs) < 2L) {
    return(NULL)
  }
  grDevices::pdf(path, width = 8, height = 6, bg = "white")
  device <- grDevices::dev.cur()
  plotted <- tryCatch(
    {
      WGCNA::plotEigengeneNetworks(
        WGCNA::orderMEs(network$MEs),
        setLabels = paste(label, "modules"),
        marDendro = c(0, 4, 1, 2),
        marHeatmap = c(3, 4, 1, 2),
        cex.lab = 0.8,
        xLabelsAngle = 90,
        plotPreservation = "standard"
      )
      TRUE
    },
    error = function(error) {
      message(
        "Could not render ", label, " eigengene network: ",
        conditionMessage(error)
      )
      FALSE
    },
    finally = {
      open_devices <- grDevices::dev.list()
      if (!is.null(open_devices) && device %in% open_devices) {
        grDevices::dev.off(device)
      }
    }
  )
  if (!plotted) {
    if (file.exists(path)) unlink(path)
    return(NULL)
  }
  normalizePath(path, mustWork = TRUE)
}

run_network_construction <- function(session_dir,
                                     method = c("blockwise", "stepwise"),
                                     max_block_size = 20000,
                                     tom_type = "auto",
                                     correlation_type = "pearson",
                                     merge_cut_height = 0.3,
                                     deep_split = 2,
                                     min_module_size = 30,
                                     random_seed = 42,
                                     kmeans_refinement = TRUE,
                                     save_adjacency = TRUE,
                                     expression_data = NULL,
                                     power_threshold = NULL,
                                     network_type = NULL,
                                     blockwise_function = build_blockwise_network,
                                     stepwise_function = build_stepwise_network,
                                     kmeans_function = apply_kmeans_refinement,
                                     create_plots = TRUE,
                                     save_artifacts = TRUE,
                                     verbose = TRUE) {
  method <- match.arg(method)
  if (is.null(expression_data) || is.null(power_threshold) || is.null(network_type)) {
    inputs <- load_network_construction_inputs(session_dir)
    expression_data <- expression_data %||network% inputs$expression_data
    power_threshold <- power_threshold %||network% inputs$power_threshold
    network_type <- network_type %||network% inputs$network_type
  }
  expression_data <- validate_network_expression(expression_data)
  method <- validate_network_parameters(
    method,
    power_threshold,
    max_block_size,
    correlation_type,
    min_module_size,
    deep_split,
    merge_cut_height,
    random_seed
  )
  tom_type <- resolve_tom_type(network_type, tom_type)
  if (!is.function(blockwise_function) || !is.function(stepwise_function) ||
      !is.function(kmeans_function)) {
    stop("Network construction hooks must be functions.", call. = FALSE)
  }

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(file.path(output_dir, "TOMfiles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message(
      "Constructing a ", network_type, " network with power ",
      power_threshold, " using the ", method, " method..."
    )
  }
  built <- if (identical(method, "blockwise")) {
    blockwise_function(
      expression_data = expression_data,
      session_dir = session_dir,
      power_threshold = power_threshold,
      network_type = network_type,
      tom_type = tom_type,
      max_block_size = max_block_size,
      correlation_type = correlation_type,
      min_module_size = min_module_size,
      deep_split = deep_split,
      merge_cut_height = merge_cut_height,
      random_seed = random_seed,
      verbose = verbose
    )
  } else {
    stepwise_function(
      expression_data = expression_data,
      session_dir = session_dir,
      power_threshold = power_threshold,
      network_type = network_type,
      tom_type = tom_type,
      correlation_type = correlation_type,
      min_module_size = min_module_size,
      deep_split = deep_split,
      merge_cut_height = merge_cut_height,
      save_adjacency = save_adjacency,
      verbose = verbose
    )
  }
  network <- built$network
  network_kmeans <- if (isTRUE(kmeans_refinement)) {
    if (verbose) message("Applying k-means module refinement...")
    kmeans_function(network, expression_data, random_seed)
  } else {
    NULL
  }

  hc_counts <- network_module_counts(network, "Hierarchical clustering")
  km_counts <- if (is.null(network_kmeans)) {
    NULL
  } else {
    network_module_counts(network_kmeans, "K-means refinement")
  }
  module_counts <- if (is.null(km_counts)) hc_counts else rbind(hc_counts, km_counts)

  parameters <- list(
    method = method,
    power_threshold = power_threshold,
    network_type = network_type,
    tom_type = tom_type,
    max_block_size = as.integer(max_block_size),
    correlation_type = correlation_type,
    merge_cut_height = merge_cut_height,
    deep_split = as.integer(deep_split),
    min_module_size = as.integer(min_module_size),
    random_seed = as.integer(random_seed),
    kmeans_refinement = isTRUE(kmeans_refinement),
    save_adjacency = isTRUE(save_adjacency)
  )
  network_result <- list(
    network = network,
    network_kmeans = network_kmeans,
    method = method,
    parameters = parameters,
    tom = built$tom,
    completed_at = Sys.time()
  )

  files <- list(
    networks = if (save_artifacts) file.path(output_dir, "3_networks.rds") else NULL,
    legacy_networks = if (save_artifacts) {
      file.path(output_dir, "3_networks.Rdata")
    } else {
      NULL
    },
    report = if (save_artifacts) file.path(output_dir, "3_network_report.rds") else NULL,
    modules_hc = if (save_artifacts) {
      file.path(output_dir, "wgcna_gene_to_module.tsv")
    } else {
      NULL
    },
    modules_kmeans = if (save_artifacts && !is.null(network_kmeans)) {
      file.path(output_dir, "wgcna_km_gene_to_module.tsv")
    } else {
      NULL
    },
    module_counts_hc = NULL,
    module_counts_kmeans = NULL,
    dendrograms_hc = character(),
    eigengenes_hc = NULL,
    eigengenes_kmeans = NULL,
    log = file.path(session_dir, "logs", "3_network_construction.log")
  )

  if (save_artifacts) {
    saveRDS(network_result, files$networks)
    net <- network
    net_km <- network_kmeans
    save(net, net_km, file = files$legacy_networks)
    utils::write.table(
      data.frame(gene = names(network$colors), module = unname(network$colors)),
      files$modules_hc,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    if (!is.null(network_kmeans)) {
      utils::write.table(
        data.frame(
          gene = names(network_kmeans$colors),
          module = unname(network_kmeans$colors)
        ),
        files$modules_kmeans,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }

  if (create_plots) {
    hc_plot <- network_module_count_plot(hc_counts, "Genes per HC module")
    files$module_counts_hc <- file.path(plots_dir, "3_genes_in_modules_HC.pdf")
    ggplot2::ggsave(files$module_counts_hc, hc_plot, width = 8, height = 6, bg = "white")
    if (!is.null(km_counts)) {
      km_plot <- network_module_count_plot(km_counts, "Genes per k-means module")
      files$module_counts_kmeans <- file.path(
        plots_dir, "3_genes_in_modules_KM.pdf"
      )
      ggplot2::ggsave(
        files$module_counts_kmeans,
        km_plot,
        width = 8,
        height = 6,
        bg = "white"
      )
    }
    files$dendrograms_hc <- save_network_dendrograms(
      network, session_dir, "HC"
    )
    files$eigengenes_hc <- save_eigengene_network_plot(
      network,
      file.path(plots_dir, "3_modules_HC_eigengenes.pdf"),
      "HC"
    )
    if (!is.null(network_kmeans)) {
      files$eigengenes_kmeans <- save_eigengene_network_plot(
        network_kmeans,
        file.path(plots_dir, "3_modules_KM_eigengenes.pdf"),
        "K-means"
      )
    }
  }

  report <- list(
    method = method,
    parameters = parameters,
    samples = nrow(expression_data),
    genes = ncol(expression_data),
    blocks = length(network$blockGenes),
    modules_hc = length(setdiff(unique(network$colors), "grey")),
    modules_kmeans = if (is.null(network_kmeans)) {
      NULL
    } else {
      length(setdiff(unique(network_kmeans$colors), "grey"))
    },
    module_counts = module_counts,
    tom = built$tom,
    files = files,
    completed_at = network_result$completed_at
  )
  class(report) <- c("wgcnamp_network_report", class(report))
  if (save_artifacts) {
    saveRDS(report, files$report)
  }
  if (verbose) {
    message(
      "Network construction complete: ", report$modules_hc,
      " non-grey HC modules across ", report$blocks, " block(s)."
    )
  }
  report
}

run_network_construction_job <- function(session_dir, job_id, ...) {
  log_dir <- file.path(session_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, "3_network_construction.log")
  started_at <- Sys.time()
  write_network_task_state(
    session_dir,
    status = "running",
    job_id = job_id,
    started_at = started_at,
    log = log_path
  )

  log_connection <- file(log_path, open = "wt")
  sink(log_connection, type = "output")
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_connection)
  }, add = TRUE)

  message("Network job ", job_id, " started at ", format(started_at))
  tryCatch({
    report <- run_network_construction(session_dir = session_dir, ...)
    write_network_task_state(
      session_dir,
      status = "complete",
      job_id = job_id,
      started_at = started_at,
      completed_at = report$completed_at,
      report = report$files$report,
      log = log_path
    )
    report
  }, error = function(error) {
    message("Network construction failed: ", conditionMessage(error))
    write_network_task_state(
      session_dir,
      status = "error",
      job_id = job_id,
      started_at = started_at,
      failed_at = Sys.time(),
      error = conditionMessage(error),
      log = log_path
    )
    stop(error)
  })
}

load_network_report <- function(session_dir) {
  path <- file.path(session_dir, "out", "3_network_report.rds")
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(error) NULL)
}

load_network_tom <- function(network_result, block = NULL) {
  manifest <- network_result$tom
  if (is.null(manifest) || length(manifest$files) == 0L) {
    stop("The network result does not contain a TOM manifest.", call. = FALSE)
  }
  if (identical(manifest$method, "stepwise")) {
    rds_file <- manifest$files[grepl("[.]rds$", manifest$files, ignore.case = TRUE)]
    if (length(rds_file) != 1L || !file.exists(rds_file)) {
      stop("The step-by-step TOM RDS file is missing.", call. = FALSE)
    }
    return(readRDS(rds_file))
  }

  files <- manifest$files
  if (is.null(block)) {
    if (length(files) != 1L) {
      stop(
        "This blockwise network has ", length(files),
        " TOM blocks. Select a block explicitly for block-specific analysis.",
        call. = FALSE
      )
    }
    block <- 1L
  }
  if (!is.numeric(block) || length(block) != 1L ||
      block < 1L || block > length(files) || block %% 1 != 0) {
    stop("Invalid TOM block index.", call. = FALSE)
  }
  file <- files[[as.integer(block)]]
  if (!file.exists(file)) {
    stop("TOM block file is missing: ", file, call. = FALSE)
  }
  environment <- new.env(parent = emptyenv())
  loaded <- load(file, envir = environment)
  matrices <- loaded[vapply(
    loaded,
    function(name) is.matrix(environment[[name]]),
    logical(1)
  )]
  if (length(matrices) != 1L) {
    stop("Could not identify one TOM matrix in block file: ", file, call. = FALSE)
  }
  environment[[matrices[[1L]]]]
}

print.wgcnamp_network_report <- function(x, ...) {
  cat("WGCNAmp network construction\n")
  cat("  Method:       ", x$method, "\n", sep = "")
  cat("  Genes:        ", x$genes, "\n", sep = "")
  cat("  Blocks:       ", x$blocks, "\n", sep = "")
  cat("  HC modules:   ", x$modules_hc, "\n", sep = "")
  if (!is.null(x$modules_kmeans)) {
    cat("  K-means modules:", x$modules_kmeans, "\n", sep = "")
  }
  invisible(x)
}

`%||network%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}
