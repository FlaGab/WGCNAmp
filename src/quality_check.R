#data quality check report

missing_QC <- function(data = data){
  
  # Validate input data
  if (!is.data.frame(data)) {
    stop("The input 'data' must be a data frame.")
  }
  
  # Quality control: check for missing values
  gsg <- WGCNA::goodSamplesGenes(t(data)) # genes as columns, samples as rows
  if (!gsg$allOK) {
    warning("Some genes or samples failed the quality control check.\n")
    print(summary(gsg$goodGenes))
    print(summary(gsg$goodSamples))
  }
  
  return(gsg)
}

HC_QC <- function(data, metadata = NULL, color_col = NULL, clust_method = "average") {
  # Validate input data
  if (!is.data.frame(data)) {
    stop("The input 'data' must be a data frame.")
  }
  if (!is.null(metadata) && !is.data.frame(metadata)) {
    stop("The input 'metadata' must be a data frame or NULL.")
  }
  if (!is.null(color_col) && !(color_col %in% colnames(metadata))) {
    stop("The specified 'color_col' is not found in the metadata.")
  }
  
  # Hierarchical clustering
  htree <- hclust(dist(t(data)), method = clust_method)
  
  # Initialize plot variable
  p <- NULL
  
  # Plot the dendrogram
  if (!is.null(metadata) && !is.null(color_col)) {
    # Check sample correspondence
    if (nrow(metadata) != ncol(data)) {
      stop("The number of samples in metadata does not match the number of samples in the data.")
    }
    
    # Ensure metadata order matches data column order
    if (!all(rownames(metadata) == colnames(data))) {
      warning("Reordering metadata to match the data.")
      metadata <- metadata[match(colnames(data), rownames(metadata)), ]
    }
    
    # Assign group colors
    group_colors <- as.factor(metadata[[color_col]])
    group_labels <- as.numeric(group_colors)
    
    # Create a colored dendrogram
    dend_col <- as.dendrogram(htree)
    dend_colored <- dendextend::color_branches(
      dend_col,
      k = length(unique(group_labels))
    )
    dend_colored <- dendextend::set(
      dend_colored,
      "labels_colors",
      group_labels[htree$order]
    )
    dend_colored <- dendextend::set(dend_colored, "labels_cex", 0.8)
    
    # Plot the dendrogram
    plot(dend_colored, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "", cex = 0.6)
    legend("topright", legend = levels(group_colors), col = 1:length(levels(group_colors)), lty = 1, cex = 0.8)
  } else {
    # Basic dendrogram without coloring
    plot(htree, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "", cex = 0.6)
  }
  
  return(invisible(htree)) # Return the hclust object for further use
}

PCA_QC <- function(data, metadata = NULL, color_col = NULL){
  
  # Validate input data
  if (!is.data.frame(data)) {
    stop("The input 'data' must be a data frame.")
  }
  if (!is.null(metadata) && !is.data.frame(metadata)) {
    stop("The input 'metadata' must be a data frame or NULL.")
  }
  if (!is.null(color_col) && !(color_col %in% colnames(metadata))) {
    stop("The specified 'color_col' is not found in the metadata.")
  }
  
  # Principal Component Analysis (PCA)
  pca <- prcomp(t(data))
  pca.dat <- as.data.frame(pca$x)
  pca.var <- pca$sdev^2 # Calculate variance
  pca.var.perc <- round(pca.var / sum(pca.var) * 100, digits = 2) # Calculate % variance explained
  
  # Enhance PCA plot with metadata-based coloring
  pca.dat$sample <- rownames(pca.dat)
  if (!is.null(metadata) && !is.null(color_col)) {
    if (!setequal(rownames(metadata), colnames(data))) {
      stop("Metadata and expression data contain different sample identifiers.")
    }
    metadata <- metadata[colnames(data), , drop = FALSE]
    pca.dat[[color_col]] <- metadata[[color_col]]
  }

  if (is.null(color_col)) {
    p <- ggplot2::ggplot(pca.dat, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(ggplot2::aes(label = sample), max.overlaps = 15)
  } else {
    p <- ggplot2::ggplot(
      pca.dat,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]])
    ) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(ggplot2::aes(label = sample), max.overlaps = 15) +
      ggplot2::scale_color_manual(
        values = grDevices::rainbow(length(unique(pca.dat[[color_col]])))
      )
  }

  p <- p +
    ggplot2::labs(
      title = "PCA Plot",
      x = paste0("PC1: ", pca.var.perc[1], "% variance"),
      y = paste0("PC2: ", pca.var.perc[2], "% variance"),
      color = color_col
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

data_QC <- function(data, metadata = NULL, color_col = NULL) {
  
  # Validate input data
  if (!is.data.frame(data)) {
    stop("The input 'data' must be a data frame.")
  }
  if (!is.null(metadata) && !is.data.frame(metadata)) {
    stop("The input 'metadata' must be a data frame or NULL.")
  }
  if (!is.null(color_col) && !(color_col %in% colnames(metadata))) {
    stop("The specified 'color_col' is not found in the metadata.")
  }
  
  list(
    missing = missing_QC(data),
    clustering = HC_QC(data, metadata = metadata, color_col = color_col),
    pca = PCA_QC(data, metadata = metadata, color_col = color_col)
  )
}
