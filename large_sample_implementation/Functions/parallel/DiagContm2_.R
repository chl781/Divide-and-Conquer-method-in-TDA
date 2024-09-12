

DiagContm2_ <- function(X,m,i,j,gap1,gap2){  
  # DiagSplit is a function to derive the diagram for the DaC.
  Diag_split=DiagSplit_parallel2_(X_split=X, i1=i, j1=j,m=m,gap1,gap2)
  #####################

  return(Diag_split)
}