#normalization script

vst_normalization <- function(data){
  
  # Create a `DESeqDataSet` object, normalize expression data
  dds <- DESeqDataSetFromMatrix(
    countData = data$filtered_data, # Our prepped data frame with counts
    colData = data$metadata, # Data frame with annotation for our samples
    design = ~1 # Here we do not need to specify a model because we are not performing a differential expression analysis
  )
  
  # Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
  dds_norm <- vst(dds)
  
  # Retrieve the normalized data from the `DESeqDataSet` and transpose it for WGCNA
  normalized_counts <- t(assay(dds_norm))
  
  return(normalized_counts)
}


