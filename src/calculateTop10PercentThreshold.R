# Function to calculate threshold for the top 10% connections in a TOM
calculateTop10PercentThreshold <- function(tom) {
  # Ensure the TOM is a numeric matrix
  if (!is.matrix(tom) || !is.numeric(tom)) {
    stop("Input must be a numeric matrix.")
  }
  
  # Flatten the TOM into a vector, excluding diagonal elements (self-connections)
  TOM_values <- as.vector(tom[lower.tri(tom)])
  
  # Remove zeros or near-zero values (optional, based on your network sparsity)
  TOM_values <- TOM_values[TOM_values > 0]
  
  # Calculate the 90th percentile (threshold for top 10% of connections)
  threshold <- quantile(TOM_values, 0.9)
  
  # Return the threshold
  return(threshold)
}

calculateEdgeWeightThreshold <- function(tom, percentile) {
  # Ensure the TOM is a numeric matrix
  if (!is.matrix(tom) || !is.numeric(tom)) {
    stop("'tom' must be a numeric matrix.")
  }
  
  if (!is.numeric(percentile) || percentile > 1 || percentile < 0) {
    stop("'percentile' must be a number between 0 and 1")
  }
  
  # Flatten the TOM into a vector, excluding diagonal elements (self-connections)
  TOM_values <- as.vector(tom[lower.tri(tom)])
  
  # Remove zeros or near-zero values (optional, based on your network sparsity)
  TOM_values <- TOM_values[TOM_values > 0]
  
  # Calculate the 90th percentile (threshold for top 10% of connections)
  threshold <- quantile(TOM_values, percentile)
  
  # Return the threshold
  return(threshold)
}