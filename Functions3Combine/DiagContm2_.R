
# This file is to parallel generate the sub-features for the data.

DiagContm2_ <- function(X,m,i,j,k,range){
  # Split data into (m-1)*(m-1) parts.
  gap1=seq(range[1,1], range[1,2], length.out = m)
  gap2=seq(range[2,1], range[2,2], length.out = m)
  gap3=seq(range[3,1], range[3,2], length.out = m)
  
  
  # DiagSplit is a function to derive the diagram for the DaC.
  Diag_split = DiagSplit_parallel_(X_split=X, i1=i, j1=j, k1=k,m=m,gap1,gap2,gap3)
  #####################
  
  return(Diag_split)
}