library(org.At.tair.db)

get_TAIRsymbols <- function(genes) {
  x <- org.At.tairSYMBOL
  mapped_tairs <- mappedkeys(x)
  xx <- as.list(x[mapped_tairs])
  annotation <- vector("character", length(genes))
  for (j in seq_along(genes)) {
    i <- genes[j]
    annotation[j] <- if (!is.null(xx[[i]][1])) xx[[i]][1] else i
  }
  
  return(annotation)
}
