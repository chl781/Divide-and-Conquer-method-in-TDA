# Generate sub-features within different sub-regions.

DiagContm2_ <- function(X,m,X_split,gap1,gap2,gap3,maxscale,maxdimension){

  
  ####### Split the data into parts with index arrays
  
  i1 = array(rep(1:(m-1),(m-1)^2),c(m-1,m-1,m-1))
  j1 = array(rep(rep(1:(m-1),(m-1)),each=m-1),c(m-1,m-1,m-1))
  k1 = array(rep(1:(m-1),each=(m-1)^2),c(m-1,m-1,m-1))
  
  # DiagSplit is a function to derive the diagram for the DaC.
  Diag_split = array(mapply(X_split=X_split, i1=i1, j1=j1, k1=k1, FUN=DiagSplit_, m=m,SIMPLIFY = F),
                      c(m-1, m-1, m-1))
  # Diag_split is the prsistent diagram in different sub-regions
  
  return(Diag_split)
}