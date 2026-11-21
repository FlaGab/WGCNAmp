# Callable implementations of pipeline steps 5 through 7.

downstream_networks <- function(network_result) {
  networks <- list(hc = network_result$network)
  if (!is.null(network_result$network_kmeans)) {
    networks$kmeans <- network_result$network_kmeans
  }
  networks
}

downstream_network_labels <- function() {
  c(hc = "Hierarchical clustering", kmeans = "K-means refinement")
}

network_gene_colors <- function(network, expression_genes) {
  colors <- network$colors
  if (is.null(colors) || length(colors) != length(expression_genes)) {
    stop("The network module assignment is incompatible with expression data.", call. = FALSE)
  }
  if (is.null(names(colors))) names(colors) <- expression_genes
  if (!setequal(names(colors), expression_genes)) {
    stop("Network modules and expression data contain different genes.", call. = FALSE)
  }
  as.character(colors[expression_genes]) |>
    stats::setNames(expression_genes)
}

downstream_control_spec <- function(session_dir) {
  inputs <- load_module_trait_inputs(session_dir)
  data <- inputs$data
  expression <- data$exprData
  if (!is.matrix(expression) && !is.data.frame(expression)) {
    stop("Prepared expression data are missing.", call. = FALSE)
  }
  expression <- as.matrix(expression)
  networks <- downstream_networks(inputs$network_result)
  labels <- downstream_network_labels()
  module_tables <- lapply(networks, function(network) {
    colors <- network_gene_colors(network, colnames(expression))
    counts <- sort(table(colors), decreasing = TRUE)
    data.frame(
      module = names(counts),
      genes = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
  colors <- lapply(networks, function(network) {
    network_gene_colors(network, colnames(expression))
  })
  continuous_traits <- if (is.null(data$traits)) character() else {
    names(data$traits)[vapply(data$traits, is.numeric, logical(1))]
  }
  metadata_columns <- if (is.null(data$metadata)) character() else {
    names(data$metadata)
  }
  list(
    samples = nrow(expression),
    genes = ncol(expression),
    networks = stats::setNames(labels[names(networks)], names(networks)),
    module_tables = module_tables,
    colors = colors,
    continuous_traits = continuous_traits,
    metadata_columns = metadata_columns,
    normalized_counts_available = !is.null(data$normalized_counts)
  )
}

load_downstream_network_parameters <- function(session_dir, network_result) {
  parameters <- network_result$parameters
  if (!is.null(parameters)) {
    return(list(
      power = parameters$power_threshold %||downstream% 6,
      network_type = parameters$network_type %||downstream% "unsigned",
      correlation_type = parameters$correlation_type %||downstream% "pearson"
    ))
  }
  threshold_path <- file.path(session_dir, "out", "2_pt.Rdata")
  if (!file.exists(threshold_path)) {
    return(list(power = 6, network_type = "unsigned", correlation_type = "pearson"))
  }
  threshold_environment <- new.env(parent = emptyenv())
  load(threshold_path, envir = threshold_environment)
  list(
    power = threshold_environment$power_threshold %||downstream% 6,
    network_type = threshold_environment$network_type %||downstream% "unsigned",
    correlation_type = threshold_environment$corType %||downstream% "pearson"
  )
}

safe_module_filename <- function(value) {
  gsub("[^A-Za-z0-9_.-]+", "_", value)
}

run_hub_analysis <- function(session_dir,
                             network_type = "hc",
                             modules,
                             min_module_membership = 0.8,
                             connectivity_quantile = 0.9,
                             trait = NULL,
                             trait_method = "pearson",
                             create_plot = TRUE,
                             save_artifacts = TRUE,
                             verbose = TRUE) {
  inputs <- load_module_trait_inputs(session_dir)
  data <- inputs$data
  expression <- as.matrix(data$exprData)
  networks <- downstream_networks(inputs$network_result)
  if (!network_type %in% names(networks)) {
    stop("The selected network is unavailable.", call. = FALSE)
  }
  if (!is.character(modules) || length(modules) == 0L) {
    stop("Select at least one module.", call. = FALSE)
  }
  if (!is.numeric(min_module_membership) || length(min_module_membership) != 1L ||
      is.na(min_module_membership) || min_module_membership < 0 ||
      min_module_membership > 1) {
    stop("Minimum module membership must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(connectivity_quantile) || length(connectivity_quantile) != 1L ||
      is.na(connectivity_quantile) || connectivity_quantile < 0 ||
      connectivity_quantile > 1) {
    stop("Connectivity quantile must be between 0 and 1.", call. = FALSE)
  }

  network <- networks[[network_type]]
  colors <- network_gene_colors(network, colnames(expression))
  available_modules <- unique(unname(colors))
  unknown_modules <- setdiff(modules, available_modules)
  if (length(unknown_modules) > 0L) {
    stop("Unknown modules: ", paste(unknown_modules, collapse = ", "), call. = FALSE)
  }
  if (is.null(rownames(network$MEs)) ||
      !setequal(rownames(network$MEs), rownames(expression))) {
    stop("Network eigengenes and expression data contain different samples.", call. = FALSE)
  }
  eigengenes <- as.matrix(network$MEs[rownames(expression), , drop = FALSE])
  gene_mm <- stats::cor(expression, eigengenes, use = "pairwise.complete.obs")
  rownames(gene_mm) <- colnames(expression)

  parameters <- load_downstream_network_parameters(session_dir, inputs$network_result)
  cor_function <- if (identical(parameters$correlation_type, "bicor")) {
    "bicor"
  } else {
    "cor"
  }
  module_tables <- lapply(modules, function(module) {
    genes <- names(colors)[colors == module]
    eigengene_column <- paste0("ME", module)
    if (!eigengene_column %in% colnames(gene_mm)) {
      stop("Module eigengene is missing for module '", module, "'.", call. = FALSE)
    }
    connectivity <- if (length(genes) == 1L) {
      stats::setNames(0, genes)
    } else {
      values <- WGCNA::softConnectivity(
        datExpr = expression[, genes, drop = FALSE],
        corFnc = cor_function,
        type = parameters$network_type,
        power = parameters$power,
        blockSize = min(1500L, length(genes)),
        verbose = 0
      )
      stats::setNames(as.numeric(values), genes)
    }
    membership <- gene_mm[genes, eigengene_column]
    connectivity_cutoff <- stats::quantile(
      connectivity,
      probs = connectivity_quantile,
      na.rm = TRUE,
      names = FALSE
    )
    table <- data.frame(
      gene = genes,
      module = module,
      module_membership = as.numeric(membership),
      abs_module_membership = abs(as.numeric(membership)),
      kWithin = as.numeric(connectivity[genes]),
      connectivity_cutoff = connectivity_cutoff,
      stringsAsFactors = FALSE
    )
    table$is_hub <- table$abs_module_membership >= min_module_membership &
      table$kWithin >= connectivity_cutoff
    table
  })
  hub_table <- do.call(rbind, module_tables)
  rownames(hub_table) <- NULL

  if (!is.null(trait) && nzchar(trait)) {
    if (is.null(data$traits) || !trait %in% names(data$traits) ||
        !is.numeric(data$traits[[trait]])) {
      stop("The selected continuous trait is unavailable.", call. = FALSE)
    }
    trait_values <- data$traits[rownames(expression), trait]
    gene_significance <- stats::cor(
      expression,
      trait_values,
      use = "pairwise.complete.obs",
      method = trait_method
    )[, 1L]
    gene_significance_p <- WGCNA::corPvalueStudent(
      matrix(gene_significance, ncol = 1L),
      nrow(expression)
    )[, 1L]
    hub_table$gene_significance <- gene_significance[hub_table$gene]
    hub_table$gene_significance_p <- gene_significance_p[hub_table$gene]
  }
  hub_table <- hub_table[order(
    hub_table$module,
    !hub_table$is_hub,
    -hub_table$abs_module_membership,
    -hub_table$kWithin
  ), , drop = FALSE]

  plot <- ggplot2::ggplot(
    hub_table,
    ggplot2::aes(x = kWithin, y = abs_module_membership, color = is_hub)
  ) +
    ggplot2::geom_point(alpha = 0.72, size = 1.8) +
    ggplot2::facet_wrap(~module, scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#94a3b8", `TRUE` = "#dc2626"),
      labels = c(`FALSE` = "Other genes", `TRUE` = "Hub candidates")
    ) +
    ggplot2::labs(
      title = "Hub-gene candidates",
      x = "Intramodular connectivity (kWithin)",
      y = "Absolute module membership",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(report = NULL, hubs = NULL, gene_membership = NULL, plot = NULL)
  if (save_artifacts) {
    files$hubs <- file.path(output_dir, "5_hub_genes.tsv")
    utils::write.table(
      hub_table, files$hubs, sep = "\t", quote = FALSE, row.names = FALSE
    )
    files$gene_membership <- file.path(
      output_dir,
      if (identical(network_type, "hc")) "geneMM_hc.csv" else "geneMM_km.csv"
    )
    utils::write.csv(
      data.frame(gene = rownames(gene_mm), gene_mm, check.names = FALSE),
      files$gene_membership,
      row.names = FALSE
    )
    for (module in modules) {
      module_hubs <- hub_table$gene[hub_table$module == module & hub_table$is_hub]
      writeLines(
        module_hubs,
        file.path(
          output_dir,
          paste0("5_", safe_module_filename(network_type), "_",
                 safe_module_filename(module), "_hub_genes.txt")
        )
      )
    }
    if (create_plot) {
      files$plot <- file.path(plots_dir, "5_hub_genes.pdf")
      ggplot2::ggsave(
        files$plot,
        plot,
        width = max(8, min(16, 4 * length(modules))),
        height = max(5, 3.5 * ceiling(length(modules) / 3)),
        bg = "white"
      )
    }
  }
  report <- list(
    parameters = list(
      network_type = network_type,
      modules = modules,
      min_module_membership = min_module_membership,
      connectivity_quantile = connectivity_quantile,
      trait = trait,
      trait_method = trait_method
    ),
    genes = nrow(hub_table),
    hubs = sum(hub_table$is_hub),
    hub_table = hub_table,
    plot = plot,
    files = files,
    completed_at = Sys.time()
  )
  class(report) <- c("wgcnamp_hub_report", class(report))
  if (save_artifacts) {
    files$report <- file.path(output_dir, "5_hub_report.rds")
    report$files$report <- files$report
    saveRDS(report, files$report)
  }
  if (verbose) {
    message("Hub analysis complete: ", report$hubs, " candidates across ",
            length(modules), " module(s).")
  }
  report
}

load_hub_report <- function(session_dir) {
  path <- file.path(session_dir, "out", "5_hub_report.rds")
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(error) NULL)
}

rbind_enrichment_tables <- function(tables) {
  if (length(tables) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(table) {
    missing <- setdiff(columns, names(table))
    for (column in missing) table[[column]] <- NA
    table[, columns, drop = FALSE]
  })
  do.call(rbind, tables)
}

parse_gene_ratio <- function(values) {
  vapply(strsplit(as.character(values), "/", fixed = TRUE), function(parts) {
    if (length(parts) != 2L) return(NA_real_)
    as.numeric(parts[[1L]]) / as.numeric(parts[[2L]])
  }, numeric(1))
}

enrichment_summary_plot <- function(results, show_categories = 15L) {
  if (nrow(results) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0, y = 0,
          label = "No enriched terms passed the selected cutoffs",
          color = "#667085", size = 4.5
        ) +
        ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
        ggplot2::theme_void()
    )
  }
  selected <- do.call(rbind, lapply(split(results, results$module), function(table) {
    table <- table[order(table$p.adjust), , drop = FALSE]
    utils::head(table, show_categories)
  }))
  selected$gene_ratio_numeric <- parse_gene_ratio(selected$GeneRatio)
  selected$Description <- factor(
    selected$Description,
    levels = rev(unique(selected$Description[order(selected$p.adjust)]))
  )
  ggplot2::ggplot(
    selected,
    ggplot2::aes(
      x = gene_ratio_numeric,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~module, scales = "free_y") +
    ggplot2::scale_color_viridis_c(direction = -1, option = "C") +
    ggplot2::labs(
      title = "GO enrichment",
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Genes"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

run_module_enrichment <- function(session_dir,
                                  network_type = "hc",
                                  modules,
                                  ontology = "BP",
                                  p_adjust_method = "BH",
                                  p_value_cutoff = 0.05,
                                  q_value_cutoff = 0.2,
                                  min_gene_set_size = 10,
                                  max_gene_set_size = 500,
                                  simplify_terms = TRUE,
                                  similarity_cutoff = 0.7,
                                  show_categories = 15,
                                  save_artifacts = TRUE,
                                  create_plot = TRUE,
                                  enrich_function = NULL,
                                  simplify_function = NULL,
                                  org_db = NULL,
                                  verbose = TRUE) {
  inputs <- load_module_trait_inputs(session_dir)
  expression_genes <- colnames(inputs$data$exprData)
  networks <- downstream_networks(inputs$network_result)
  if (!network_type %in% names(networks)) {
    stop("The selected network is unavailable.", call. = FALSE)
  }
  colors <- network_gene_colors(networks[[network_type]], expression_genes)
  if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
  unknown <- setdiff(modules, unique(unname(colors)))
  if (length(unknown) > 0L) {
    stop("Unknown modules: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (!ontology %in% c("BP", "MF", "CC", "ALL")) {
    stop("Ontology must be BP, MF, CC, or ALL.", call. = FALSE)
  }
  if (is.null(enrich_function)) {
    required <- c("clusterProfiler", "org.At.tair.db")
    available <- vapply(required, function(package) {
      suppressPackageStartupMessages(requireNamespace(package, quietly = TRUE))
    }, logical(1))
    missing <- required[!available]
    if (length(missing) > 0L) {
      stop("Install enrichment packages: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    enrich_function <- clusterProfiler::enrichGO
  }
  if (is.null(org_db)) {
    org_db <- getExportedValue("org.At.tair.db", "org.At.tair.db")
  }
  if (isTRUE(simplify_terms) && is.null(simplify_function) &&
      suppressPackageStartupMessages(requireNamespace("clusterProfiler", quietly = TRUE))) {
    simplify_function <- clusterProfiler::simplify
  }

  universe <- names(colors)
  results_by_module <- list()
  for (module in modules) {
    genes <- names(colors)[colors == module]
    enrichment <- enrich_function(
      gene = genes,
      OrgDb = org_db,
      keyType = "TAIR",
      ont = ontology,
      pvalueCutoff = p_value_cutoff,
      pAdjustMethod = p_adjust_method,
      universe = universe,
      qvalueCutoff = q_value_cutoff,
      minGSSize = as.integer(min_gene_set_size),
      maxGSSize = as.integer(max_gene_set_size),
      readable = FALSE
    )
    table <- as.data.frame(enrichment)
    if (isTRUE(simplify_terms) && nrow(table) > 0L &&
        !is.null(simplify_function) && !identical(ontology, "ALL")) {
      enrichment <- simplify_function(
        enrichment,
        cutoff = similarity_cutoff,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang"
      )
      table <- as.data.frame(enrichment)
    }
    if (nrow(table) > 0L) {
      table$module <- module
      results_by_module[[module]] <- table
    }
  }
  combined <- rbind_enrichment_tables(results_by_module)
  rownames(combined) <- NULL
  plot <- enrichment_summary_plot(combined, as.integer(show_categories))

  output_dir <- file.path(session_dir, "out", "GO")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(report = NULL, combined = NULL, module_tables = list(), plot = NULL)
  if (save_artifacts) {
    files$combined <- file.path(output_dir, "6_go_enrichment.tsv")
    utils::write.table(
      combined, files$combined, sep = "\t", quote = FALSE, row.names = FALSE
    )
    for (module in names(results_by_module)) {
      path <- file.path(
        output_dir,
        paste0("6_", safe_module_filename(module), "_", ontology, "_enrichment.csv")
      )
      utils::write.csv(results_by_module[[module]], path, row.names = FALSE)
      files$module_tables[[module]] <- path
    }
    if (create_plot) {
      files$plot <- file.path(plots_dir, "6_go_enrichment.pdf")
      ggplot2::ggsave(
        files$plot,
        plot,
        width = 11,
        height = max(6, min(20, 4 * ceiling(max(1, length(modules)) / 2))),
        bg = "white"
      )
    }
  }
  report <- list(
    parameters = list(
      network_type = network_type,
      modules = modules,
      ontology = ontology,
      p_adjust_method = p_adjust_method,
      p_value_cutoff = p_value_cutoff,
      q_value_cutoff = q_value_cutoff,
      min_gene_set_size = min_gene_set_size,
      max_gene_set_size = max_gene_set_size,
      simplify_terms = simplify_terms,
      similarity_cutoff = similarity_cutoff,
      show_categories = show_categories
    ),
    terms = nrow(combined),
    modules_with_terms = length(results_by_module),
    results = combined,
    plot = plot,
    files = files,
    completed_at = Sys.time()
  )
  class(report) <- c("wgcnamp_enrichment_report", class(report))
  if (save_artifacts) {
    files$report <- file.path(output_dir, "6_enrichment_report.rds")
    report$files$report <- files$report
    saveRDS(report, files$report)
  }
  if (verbose) {
    message("Enrichment complete: ", report$terms, " terms across ",
            report$modules_with_terms, " module(s).")
  }
  report
}

load_enrichment_report <- function(session_dir) {
  path <- file.path(session_dir, "out", "GO", "6_enrichment_report.rds")
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(error) NULL)
}

scale_expression_rows <- function(matrix) {
  scaled <- t(scale(t(matrix)))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

expression_long_table <- function(matrix) {
  data.frame(
    gene = rep(rownames(matrix), each = ncol(matrix)),
    sample = rep(colnames(matrix), times = nrow(matrix)),
    value = as.vector(t(matrix)),
    stringsAsFactors = FALSE
  )
}

run_gene_visualization <- function(session_dir,
                                   network_type = "hc",
                                   modules,
                                   genes = character(),
                                   plot_type = c("gene_profiles", "module_summary", "eigengene_heatmap"),
                                   group_column = NULL,
                                   scale_rows = TRUE,
                                   create_plot = TRUE,
                                   save_artifacts = TRUE,
                                   verbose = TRUE) {
  plot_type <- match.arg(plot_type)
  inputs <- load_module_trait_inputs(session_dir)
  data <- inputs$data
  networks <- downstream_networks(inputs$network_result)
  if (!network_type %in% names(networks)) {
    stop("The selected network is unavailable.", call. = FALSE)
  }
  if (is.null(data$normalized_counts)) {
    stop("Normalized counts are unavailable. Re-run data preparation.", call. = FALSE)
  }
  counts <- as.matrix(data$normalized_counts)
  metadata <- data$metadata
  if (!setequal(colnames(counts), rownames(metadata))) {
    stop("Normalized counts and metadata contain different samples.", call. = FALSE)
  }
  metadata <- metadata[colnames(counts), , drop = FALSE]
  network <- networks[[network_type]]
  colors <- network_gene_colors(network, rownames(counts))
  if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
  unknown_modules <- setdiff(modules, unique(unname(colors)))
  if (length(unknown_modules) > 0L) {
    stop("Unknown modules: ", paste(unknown_modules, collapse = ", "), call. = FALSE)
  }
  module_genes <- names(colors)[colors %in% modules]
  if (length(genes) == 0L) genes <- utils::head(module_genes, 12L)
  unknown_genes <- setdiff(genes, module_genes)
  if (length(unknown_genes) > 0L) {
    stop("Selected genes are not in the selected modules: ",
         paste(utils::head(unknown_genes, 10L), collapse = ", "), call. = FALSE)
  }
  if (length(genes) > 100L) {
    stop("Select at most 100 genes for one visualization.", call. = FALSE)
  }
  if (!is.null(group_column) && nzchar(group_column) &&
      !group_column %in% names(metadata)) {
    stop("The selected grouping column is unavailable.", call. = FALSE)
  }
  group_column <- if (is.null(group_column) || !nzchar(group_column)) NULL else group_column
  sample_order <- rownames(metadata)

  if (identical(plot_type, "eigengene_heatmap")) {
    eigengene_columns <- intersect(paste0("ME", modules), colnames(network$MEs))
    values <- t(as.matrix(network$MEs[sample_order, eigengene_columns, drop = FALSE]))
    result_table <- data.frame(
      module = rep(sub("^ME", "", rownames(values)), each = ncol(values)),
      sample = rep(colnames(values), times = nrow(values)),
      value = as.vector(t(values)),
      stringsAsFactors = FALSE
    )
    result_table$sample <- factor(result_table$sample, levels = sample_order)
    plot <- ggplot2::ggplot(
      result_table,
      ggplot2::aes(x = sample, y = module, fill = value)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::scale_fill_gradient2(
        low = "#2563eb", mid = "white", high = "#dc2626", midpoint = 0,
        name = "Eigengene"
      ) +
      ggplot2::labs(title = "Module eigengenes", x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)
      )
  } else if (identical(plot_type, "module_summary")) {
    module_matrices <- lapply(modules, function(module) {
      selected <- names(colors)[colors == module]
      matrix <- counts[selected, sample_order, drop = FALSE]
      if (isTRUE(scale_rows)) matrix <- scale_expression_rows(matrix)
      colMeans(matrix, na.rm = TRUE)
    })
    names(module_matrices) <- modules
    result_table <- do.call(rbind, lapply(modules, function(module) {
      data.frame(
        module = module,
        sample = sample_order,
        value = as.numeric(module_matrices[[module]]),
        stringsAsFactors = FALSE
      )
    }))
    result_table$sample <- factor(result_table$sample, levels = sample_order)
    if (!is.null(group_column)) {
      result_table$group <- as.character(metadata[as.character(result_table$sample), group_column])
    } else {
      result_table$group <- result_table$module
    }
    plot <- ggplot2::ggplot(
      result_table,
      ggplot2::aes(x = sample, y = value, group = module, color = module)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::facet_wrap(~module, scales = "free_y") +
      ggplot2::labs(
        title = "Module expression summaries", x = NULL,
        y = if (isTRUE(scale_rows)) "Mean gene z-score" else "Mean normalized count",
        color = "Module"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
        legend.position = "none"
      )
  } else {
    matrix <- counts[genes, sample_order, drop = FALSE]
    if (isTRUE(scale_rows)) matrix <- scale_expression_rows(matrix)
    result_table <- expression_long_table(matrix)
    result_table$module <- unname(colors[result_table$gene])
    result_table$sample <- factor(result_table$sample, levels = sample_order)
    if (!is.null(group_column)) {
      result_table$group <- as.character(metadata[as.character(result_table$sample), group_column])
    } else {
      result_table$group <- "Samples"
    }
    plot <- ggplot2::ggplot(
      result_table,
      ggplot2::aes(x = sample, y = value, group = 1, color = group)
    ) +
      ggplot2::geom_line(linewidth = 0.65) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::facet_wrap(~gene, scales = "free_y") +
      ggplot2::labs(
        title = "Selected gene expression", x = NULL,
        y = if (isTRUE(scale_rows)) "Gene z-score" else "Normalized count",
        color = group_column %||downstream% "Group"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  }

  output_dir <- file.path(session_dir, "out")
  plots_dir <- file.path(session_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(report = NULL, data = NULL, plot = NULL)
  if (save_artifacts) {
    files$data <- file.path(output_dir, "7_visualization_data.tsv")
    utils::write.table(
      result_table, files$data, sep = "\t", quote = FALSE, row.names = FALSE
    )
    if (create_plot) {
      files$plot <- file.path(plots_dir, "7_gene_visualization.pdf")
      ggplot2::ggsave(
        files$plot,
        plot,
        width = 11,
        height = if (identical(plot_type, "gene_profiles")) {
          max(6, 2.8 * ceiling(length(genes) / 3))
        } else {
          max(6, 3 * ceiling(length(modules) / 3))
        },
        bg = "white",
        limitsize = FALSE
      )
    }
  }
  report <- list(
    parameters = list(
      network_type = network_type,
      modules = modules,
      genes = genes,
      plot_type = plot_type,
      group_column = group_column,
      scale_rows = isTRUE(scale_rows)
    ),
    rows = nrow(result_table),
    plot = plot,
    data = result_table,
    files = files,
    completed_at = Sys.time()
  )
  class(report) <- c("wgcnamp_visualization_report", class(report))
  if (save_artifacts) {
    files$report <- file.path(output_dir, "7_visualization_report.rds")
    report$files$report <- files$report
    saveRDS(report, files$report)
  }
  if (verbose) message("Gene visualization complete.")
  report
}

load_visualization_report <- function(session_dir) {
  path <- file.path(session_dir, "out", "7_visualization_report.rds")
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(error) NULL)
}

`%||downstream%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}
