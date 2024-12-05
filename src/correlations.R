#this script create binary contrast matrix for regression of categorical data

function binary_contrast_matrix(data = data$metadata,
                                variables){
  if (all(variables %in% data)){
    
    bin.mat <- binarizeCategoricalColumns(data,
                                          #convertColumns = ,
                                          considerColumns = variables,
                                          #maxOrdinalLevels = ,
                                          #levelOrder = ,
                                          minCount = 1,
                                          includePairwise = TRUE,
                                          includeLevelVsAll = FALSE,
                                          dropUninformative = TRUE,
                                          includePrefix = FALSE,
                                          #prefixSep = " vs. ",
                                          levelSep = " vs. "
    )
    
    print(paste("Binary matrix includes:", colnames(bin.mat), collapse = "/n"))
  } else {
    stop("indicate variables that are included in the data provided")
  }
}