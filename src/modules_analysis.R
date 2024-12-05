#functions
hub_genes <- function(network, module, exprData = NULL, TOM = TOM, q_kWithin = 0.9, minMM = 0.8){
  
  if(!(module %in% network$colors)){
    stop("The specified module does not exist in the provided network. Please check the module name.")
  } else {
    if(!exists("geneMM", where = network)){
      if(is.null(exprData)){
        stop("Gene module membership matrix not found. Please provide 'exprData' to calculate it.")
        } else {
          warning("Gene module membership matrix not found. Calculating it...")
          network$geneMM <- cor(exprData, network$MEs)
        }
      }
    if(!exists("IMconnectivity", where = network)){
      if(is.null(TOM)){
        stop("Intramodular connectivity not found. Please provide 'TOM' to calculate it.")
        } else {
          warning("Intramodular connectivity not found. Calculating it...")
          network$IMconnectivity <- intramodularConnectivity(adjMat = TOM, colors = network$colors)
        }
      }
    
    module_hubs <- network$IMconnectivity[
      network$IMconnectivity$kWithin > quantile(network$IMconnectivity$kWithin, q_kWithin) &
      abs(network$geneMM[, paste0("ME", module)]) > minMM &
      rownames(network$IMconnectivity) %in% names(network$colors[network$colors == module]),
    ]
  }
  
  return(module_hubs)
}
