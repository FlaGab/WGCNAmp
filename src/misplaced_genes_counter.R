#misplaced genes counter function
misplaced.genes.dataframe <- function(MM = genModuleMembership, P = bwnet_km$colors){
  correct = 0
  misplaced = 0
  for(g in rownames(MM)){
    if (max(MM[g,]) == MM[g, paste0("ME", P[g])]){
      #print(paste0(g, " : correctly assigned module"))
      correct = correct + 1
    }
    else{
      #print(paste0(g, " : misplaced")) 
      misplaced = misplaced +1
    }
  }
  out = c(correct = correct, misplaced = misplaced)
  
  return(out)
}
