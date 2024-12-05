#convert colors vector to numerical vedctor

col2num <- function(color_vector){
  
  numcode <- 0:c(length(unique(color_vector))-1)
  names(numcode) <- unique(color_vector)
  result <- c()
  for (i in color_vector){
    result <- c(result, numcode[i])
  }
  names(result) <- names(color_vector)
  
  return(result)
}