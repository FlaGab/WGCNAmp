#data filtering functions

#function for excluding outlier samples

exclude_samples <- function(data, samples = NULL){
  
  if(!is.null(samples)){
    if(!is.vector(samples)){
      stop("The input 'samples' must be a vector. select a subset of:\n ", paste(rownames(data$metadata), collapse = "\n"))
    }
    if(!(samples %in% rownames(data$metadata))){
      stop("you must specify samples in the data submitted. select a subset of:\n ", paste(rownames(data$metadata), collapse = "\n"))
    }
    else{
      out <- list (counts = data$data_no_missing[, !(colnames(data$data_no_missing) %in% samples)],
                   metadata = data$metadata[!(rownames(data$metadata) %in% samples), ],
                   traits = filter(data$traits, rownames(data$traits) != samples),#data$traits[!(rownames(data$traits) %in% samples), ],
                   gene_length = data$gene_length)
    }
  }
  
  return(out)
}

#function for filtering low expressed genes 

#option 1 - filtering genes with less than n total reads across all samples

#option 2 - filtering genes with less than n reads in m saples

#option 3 - filtering genes with less than n reads in p percent of the samples

low_reads_filtering <- function(data, method = c("1", "2", "3"), n = 0 , m = 0, p = 0){
  
  if (!is.data.frame(data)) {
    stop("The input 'data' must be a data frame containing raw reads count.")
  }
  
  if (!(method %in% c("1", "2", "3"))){
    
    stop("Invalid method. Please choose one of the following:\n
          1 - filtering genes with less than n total reads across all samples;\n
          2 - filtering genes with less than n reads in m saples;\n
          3 - filtering genes with less than n reads in p percent of the samples.")
  }
  
  # Validate n
  if (!is.numeric(n) || n < 0) {
    stop("'n' must be a non-negative numeric value.")
  }
  
  if (method == "1"){
    print(paste0("Filtering by method ",method,"."))
    print(paste0("Removing genes with less than ", n ," total reads across all samples"))
    filtered_data <- data[rowSums(data) >= n, ]
  }
  
  else if (method == "2"){
    
    # Validate m
    if (!is.numeric(m) || m <= 0 || m > ncol(data)) {
      stop("'m' must be a numeric value greater than 0 and less than or equal to the number of samples.")
    }
    
    print(paste0("Filtering by method ",method,"."))
    print(paste0("Removing genes with less than ", n ," reads in ", m, " or more samples"))
    filtered_data <- data[rowSums(data > n) >= m, ]
  }
  
  else if (method == "3") {
    # Validate p
    if (!is.numeric(p) || p <= 0 || p > 1) {
      stop("'p' must be a numeric value between 0 and 1 (proportion of samples).")
    }
    
    print(paste0("Filtering by method ",method,"."))
    print(paste0("Removing genes with less than ", n ," reads in ", p*100, "% of the samples or more"))
    filtered_data <- data[rowSums(data > n) >= ncol(data)*p,]
  }
  
  removed_genes <- nrow(data)-nrow(filtered_data)
  perc_removed_genes <- round((removed_genes/nrow(data))*100, digits = 2)
  print(paste0("Total removed genes: ", removed_genes, " (",perc_removed_genes,"%)"))
  return(filtered_data)
}