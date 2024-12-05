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

binary_contrast_matrix <- function(data, variables, vsAll = TRUE, verbose = TRUE) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (is.null(variables) || length(variables) == 0) {
    stop("Please specify at least one variable for binarization.")
  }
  if (!all(variables %in% colnames(data))) {
    missing_vars <- variables[!variables %in% colnames(data)]
    stop(paste("The following variables are not in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Generate binary matrix
  bin.mat <- binarizeCategoricalColumns(data,
                                        considerColumns = variables,
                                        minCount = 1,
                                        includePairwise = TRUE,
                                        includeLevelVsAll = vsAll,
                                        dropUninformative = TRUE,
                                        includePrefix = TRUE,
                                        prefixSep = ": ",
                                        levelSep = " vs. ")
  rownames(bin.mat) <- rownames(data)
  
  # Print summary if verbose is TRUE
  if (verbose) {
    message("Binary matrix generated with dimensions: ", nrow(bin.mat), " rows and ", ncol(bin.mat), " columns.")
    message("Contrasts:")
    cat(paste(seq_along(colnames(bin.mat)), colnames(bin.mat), sep = ". "), sep = "\n")
  }
  
  # Return the binary matrix
  return(bin.mat)
}

correlation_catgorical_var <- function(network = net, binary_matrix = binary.matrix){
  
  #define number of genes and samples
  nSamples = nrow(binary_matrix)
  nGenes = length(network$colors)
  
  module.trait.corr <- cor(network$MEs, 
                           binary_matrix,
                           use = "p",
                           method = "p") #change the method based on your data
  module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
  
  #visualize module-trait association
  heatmap.data <- merge(network$MEs, binary_matrix, by = "row.names")
  rownames(heatmap.data) <- heatmap.data$Row.names
  heatmap.data$Row.names <- NULL
  
  CLP <- CorLevelPlot(heatmap.data,
                      x = names(heatmap.data)[(ncol(heatmap.data)-ncol(binary.matrix)+1):ncol(heatmap.data)],
                      y = names(heatmap.data)[1:(ncol(heatmap.data)-ncol(binary.matrix))],
                      col = c("blue1", "skyblue", "white", "pink", "red"),
                      rotLabX = 30,
                      plotRsquared = F
  )
  
  return(list(r_vals = module.trait.corr, p_vals = module.trait.corr.pvals, hm = CLP))
}

correlation_categorical_var_lm <- function(network = net, 
                                           binary_matrix = binary.matrix, 
                                           plot_type = "dotplot") {
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
