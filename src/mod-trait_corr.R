#module-reai
generate_random_df <- function(n, m, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed) # Ensure reproducibility if a seed is provided
  }
  
  # Generate n rows and m columns of random values
  random_matrix <- matrix(runif(n * m), nrow = n, ncol = m)
  
  # Normalize rows to sum to 1
  random_matrix <- t(apply(random_matrix, 1, function(row) row / sum(row)))
  
  # Convert to a data frame
  random_df <- as.data.frame(random_matrix)
  
  # Optionally name the columns
  colnames(random_df) <- paste0("Var", seq_len(m))
  
  return(random_df)
}

continuous_traits_correlations <- function(data = data, network = net, method = "pearson"){
  #define number of genes and samples
  nSamples = nrow(data$traits)
  nGenes = ncol(data$exprData)
  
  module.trait.corr <- cor(network$MEs, data$traits,
                           method = method)
  module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
  
  #visualize module-trait association
  heatmap.data <- merge(network$MEs, data$traits, by = "row.names")
  rownames(heatmap.data) <- heatmap.data$Row.names
  heatmap.data$Row.names <- NULL
  
  CLP <- CorLevelPlot(heatmap.data,
                      x = names(heatmap.data)[(length(unique(network$colors))+1):(length(unique(network$colors))+length(colnames(data$traits)))],
                      y = names(heatmap.data)[1:length(unique(network$colors))],
                      col = c("blue1", "skyblue", "white", "pink", "red"),
                      rotLabX = 45,
                      cexCorval = 0.7
  )
  
  return(CLP)
}

detect_categorical_variables <- function(data) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  names(data)[vapply(
    data,
    function(column) {
      is.factor(column) || is.character(column) || is.logical(column)
    },
    logical(1)
  )]
}

prepare_categorical_metadata <- function(data,
                                         variables = NULL,
                                         reference_levels = NULL) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  if (is.null(variables)) {
    variables <- detect_categorical_variables(data)
  }

  if (!is.character(variables) || length(variables) == 0L) {
    stop(
      "No categorical variables were selected. Supply column names in ",
      "'variables', or convert categorical metadata columns to factors."
    )
  }

  variables <- unique(variables)

  if (!all(variables %in% colnames(data))) {
    missing_vars <- variables[!variables %in% colnames(data)]
    stop(
      "The following categorical variables are not in the data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (!is.null(reference_levels)) {
    reference_levels <- unlist(reference_levels, use.names = TRUE)

    if (is.null(names(reference_levels)) ||
        any(names(reference_levels) == "") ||
        anyDuplicated(names(reference_levels))) {
      stop("'reference_levels' must be a uniquely named vector or list.")
    }

    unknown_reference_variables <- setdiff(names(reference_levels), variables)
    if (length(unknown_reference_variables) > 0L) {
      stop(
        "Reference levels were supplied for unselected variables: ",
        paste(unknown_reference_variables, collapse = ", ")
      )
    }
  }

  prepared_data <- data
  observed_level_counts <- vapply(
    variables,
    function(variable) {
      nlevels(droplevels(factor(prepared_data[[variable]])))
    },
    integer(1)
  )
  skipped_variables <- variables[observed_level_counts < 2L]
  variables <- variables[observed_level_counts >= 2L]

  if (length(variables) == 0L) {
    stop(
      "No informative categorical variables remain: every selected variable ",
      "has fewer than two observed levels."
    )
  }

  for (variable in variables) {
    values <- droplevels(factor(prepared_data[[variable]]))

    if (!is.null(reference_levels) && variable %in% names(reference_levels)) {
      reference <- unname(reference_levels[[variable]])

      if (length(reference) != 1L || is.na(reference) ||
          !(reference %in% levels(values))) {
        stop(
          "Reference level '", reference, "' is not an observed level of '",
          variable, "'. Available levels: ",
          paste(levels(values), collapse = ", ")
        )
      }

      values <- stats::relevel(values, ref = reference)
    }

    prepared_data[[variable]] <- values
  }

  list(
    data = prepared_data,
    variables = variables,
    skipped_variables = skipped_variables
  )
}

drop_redundant_binary_contrasts <- function(binary_matrix, variables) {
  if (ncol(binary_matrix) < 2L) {
    return(binary_matrix)
  }

  contrast_variables <- vapply(
    colnames(binary_matrix),
    function(contrast) {
      matches <- variables[startsWith(contrast, paste0(variables, ": "))]
      if (length(matches) == 0L) NA_character_ else matches[[1L]]
    },
    character(1)
  )

  redundant <- rep(FALSE, ncol(binary_matrix))

  for (variable in unique(stats::na.omit(contrast_variables))) {
    positions <- which(contrast_variables == variable)
    signatures <- vapply(
      positions,
      function(position) {
        paste(binary_matrix[, position], collapse = "\r")
      },
      character(1)
    )
    redundant[positions] <- duplicated(signatures)
  }

  dropped <- colnames(binary_matrix)[redundant]
  result <- binary_matrix[, !redundant, drop = FALSE]
  attr(result, "dropped_redundant_contrasts") <- dropped
  result
}

binary_contrast_matrix <- function(data,
                                   variables = NULL,
                                   reference_levels = NULL,
                                   vsAll = TRUE,
                                   includePairwise = TRUE,
                                   minCount = 1,
                                   drop_redundant = TRUE,
                                   verbose = TRUE) {
  prepared <- prepare_categorical_metadata(
    data = data,
    variables = variables,
    reference_levels = reference_levels
  )

  if (!is.logical(vsAll) || length(vsAll) != 1L || is.na(vsAll)) {
    stop("'vsAll' must be TRUE or FALSE.")
  }
  if (!is.logical(includePairwise) ||
      length(includePairwise) != 1L ||
      is.na(includePairwise)) {
    stop("'includePairwise' must be TRUE or FALSE.")
  }
  if (!vsAll && !includePairwise) {
    stop("At least one contrast type must be enabled.")
  }

  # Generate binary matrix
  bin.mat <- WGCNA::binarizeCategoricalColumns(
    prepared$data[, prepared$variables, drop = FALSE],
    convertColumns = prepared$variables,
    considerColumns = prepared$variables,
    minCount = minCount,
    includePairwise = includePairwise,
    includeLevelVsAll = vsAll,
    dropUninformative = TRUE,
    includePrefix = TRUE,
    prefixSep = ": ",
    levelSep = " vs. "
  )
  rownames(bin.mat) <- rownames(data)

  if (ncol(bin.mat) == 0L) {
    stop("No informative binary contrasts could be generated.")
  }

  if (drop_redundant) {
    bin.mat <- drop_redundant_binary_contrasts(
      binary_matrix = bin.mat,
      variables = prepared$variables
    )
  }
  attr(bin.mat, "skipped_single_level_variables") <- prepared$skipped_variables

  # Print summary if verbose is TRUE
  if (verbose) {
    if (length(prepared$skipped_variables) > 0L) {
      message(
        "Skipped categorical variables with fewer than two observed levels: ",
        paste(prepared$skipped_variables, collapse = ", ")
      )
    }

    message("Categorical variables: ", paste(prepared$variables, collapse = ", "))
    for (variable in prepared$variables) {
      message(
        "Reference level for '", variable, "': ",
        levels(prepared$data[[variable]])[[1L]]
      )
    }

    dropped <- attr(bin.mat, "dropped_redundant_contrasts")
    if (length(dropped) > 0L) {
      message(
        "Dropped redundant within-variable contrasts: ",
        paste(dropped, collapse = ", ")
      )
    }

    message(
      "Binary matrix generated with dimensions: ",
      nrow(bin.mat), " rows and ", ncol(bin.mat), " columns."
    )
    message("Contrasts:")
    cat(paste0("- ", colnames(bin.mat)), sep = "\n")
  }

  # Return the binary matrix
  return(bin.mat)
}

select_binary_contrasts <- function(binary_matrix,
                                    include = NULL,
                                    exclude = NULL,
                                    verbose = TRUE) {
  if (!is.matrix(binary_matrix) && !is.data.frame(binary_matrix)) {
    stop("'binary_matrix' must be a matrix or data frame.")
  }

  available <- colnames(binary_matrix)
  if (is.null(available)) {
    stop("'binary_matrix' must have contrast names.")
  }

  if (!is.null(include)) {
    missing_includes <- setdiff(include, available)
    if (length(missing_includes) > 0L) {
      stop(
        "Requested contrasts were not generated: ",
        paste(missing_includes, collapse = ", "),
        ". Available contrasts: ", paste(available, collapse = ", ")
      )
    }
    selected <- unique(include)
  } else {
    selected <- available
  }

  if (!is.null(exclude)) {
    unknown_excludes <- setdiff(exclude, available)
    if (length(unknown_excludes) > 0L) {
      stop(
        "Contrasts requested for exclusion were not generated: ",
        paste(unknown_excludes, collapse = ", ")
      )
    }
    selected <- setdiff(selected, exclude)
  }

  if (length(selected) == 0L) {
    stop("Contrast selection removed every generated contrast.")
  }

  skipped_variables <- attr(
    binary_matrix,
    "skipped_single_level_variables"
  )
  result <- binary_matrix[, selected, drop = FALSE]
  attr(result, "skipped_single_level_variables") <- skipped_variables

  if (verbose) {
    message("Selected contrasts: ", paste(selected, collapse = ", "))
  }

  result
}

align_binary_matrix <- function(network, binary_matrix) {
  eigengene_samples <- rownames(network$MEs)
  contrast_samples <- rownames(binary_matrix)

  if (is.null(eigengene_samples) || is.null(contrast_samples)) {
    stop("Network eigengenes and the binary matrix must have sample row names.")
  }

  if (!setequal(eigengene_samples, contrast_samples)) {
    missing_from_contrasts <- setdiff(eigengene_samples, contrast_samples)
    missing_from_network <- setdiff(contrast_samples, eigengene_samples)
    stop(
      "Network eigengenes and binary contrasts contain different samples. ",
      "Missing from contrasts: ",
      paste(missing_from_contrasts, collapse = ", "),
      "; missing from network: ",
      paste(missing_from_network, collapse = ", ")
    )
  }

  binary_matrix[eigengene_samples, , drop = FALSE]
}

correlation_categorical_var <- function(network, binary_matrix){
  binary_matrix <- align_binary_matrix(network, binary_matrix)
  
  #define number of genes and samples
  nSamples = nrow(binary_matrix)
  nGenes = length(network$colors)
  
  module.trait.corr <- cor(network$MEs, 
                           binary_matrix,
                           use = "p",
                           method = "pearson")
  module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
  
  #visualize module-trait association
  heatmap.data <- merge(network$MEs, binary_matrix, by = "row.names")
  rownames(heatmap.data) <- heatmap.data$Row.names
  heatmap.data$Row.names <- NULL
  
  CLP <- CorLevelPlot(heatmap.data,
                      x = names(heatmap.data)[(ncol(heatmap.data)-ncol(binary_matrix)+1):ncol(heatmap.data)],
                      y = names(heatmap.data)[1:(ncol(heatmap.data)-ncol(binary_matrix))],
                      col = c("blue1", "skyblue", "white", "pink", "red"),
                      rotLabX = 30,
                      plotRsquared = F
  )
  
  return(list(r_vals = module.trait.corr, p_vals = module.trait.corr.pvals, hm = CLP))
}

# Backward-compatible alias for the original misspelled function name.
correlation_catgorical_var <- correlation_categorical_var

correlation_categorical_var_lm <- function(network,
                                           binary_matrix,
                                           plot_type = "dotplot") {
  binary_matrix <- align_binary_matrix(network, binary_matrix)

  # Initialize an empty list to store results
  results <- list()
  
  # Loop over contrasts (columns of binary_matrix)
  for (i in 1:ncol(binary_matrix)) {
    contrast_name <- colnames(binary_matrix)[i] # Get contrast name
    
    # Initialize a list to store module results for this contrast
    contrast_results <- list()
    
    # Loop over modules (columns of network$MEs)
    for (j in 1:ncol(network$MEs)) {
      module_name <- colnames(network$MEs)[j] # Get module name
      eigengene <- network$MEs[, j]          # Extract eigengene values
      
      # Fit linear model
      lm_model <- lm(eigengene ~ binary_matrix[[i]])
      
      # Store important metrics from the summary
      summary_model <- summary(lm_model)
      contrast_results[[module_name]] <- list(
        estimate = coef(summary_model)[2, "Estimate"],
        p_value = coef(summary_model)[2, "Pr(>|t|)"],
        r_squared = summary_model$r.squared,
        adj_r_squared = summary_model$adj.r.squared
      )
    }
    
    # Store the results for this contrast in the main results list
    results[[contrast_name]] <- contrast_results
  }
  
  # Combine results into a data frame for plotting
  plot_data <- do.call(rbind, lapply(names(results), function(contrast) {
    do.call(rbind, lapply(names(results[[contrast]]), function(module) {
      data.frame(
        Contrast = contrast,
        Module = module,
        Estimate = results[[contrast]][[module]]$estimate,
        P_Value = results[[contrast]][[module]]$p_value
      )
    }))
  }))
  
  # Add a significance column with asterisks
  plot_data$Significance <- with(plot_data, ifelse(
    P_Value < 0.001, "***",
    ifelse(P_Value < 0.01, "**",
           ifelse(P_Value < 0.05, "*", ""))
  ))
    
  # Create plots based on plot_type
  if (plot_type == "heatmap") {
    p <- ggplot(plot_data, aes(x = Contrast, y = Module)) +
      geom_tile(aes(fill = Estimate), color = "white") +
      geom_text(aes(label = paste0(round(Estimate, 2), "\n", Significance)), size = 4) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      labs(x = "Contrast", y = "Module", fill = "Mean Eigengenes difference") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else if (plot_type == "dotplot") {
    p <- ggplot(plot_data, aes(x = Contrast, y = Module)) +
      geom_point(aes(size = -log(P_Value), color = Estimate)) +
      geom_text(aes(label = Significance), vjust = -1) +
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      scale_size_continuous(range = c(3, 10)) +
      labs(x = "Contrast", y = "Module", size = "-log(p-value)", color = "Mean Eigengenes difference") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    stop("Invalid plot_type. Use 'heatmap' or 'dotplot'.")
  }
  
  # Print the plot
  print(p)
  
  # Return the results as a structured object
  return(results)
}
