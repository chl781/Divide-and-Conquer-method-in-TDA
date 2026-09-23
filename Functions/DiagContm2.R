

DiagContm2 <- function(X,m,maxscale,maxdimension,range){
  # Split data into (m-1)*(m-1) parts.
  gap1=seq(range[1,1], range[1,2], length.out = m)
  gap2=seq(range[2,1], range[2,2], length.out = m)
  
  #
  X_split = matrix(list(),m-1,m-1)
  for (i in 1:(m-1)) {
    for(j in 1:(m-1)){
      X_split[[i,j]] = X[ X[,1] >= gap1[i] & X[,2] >= gap2[j] 
                         & X[,1] < gap1[i+1] & X[,2] < gap2[j+1],]
    }
  }
  
  ##################### Begin to split the data into parts.
  
  i1 = matrix(rep(1:(m-1),m-1),m-1,m-1)
  j1 = matrix(rep(1:(m-1),m-1),m-1,m-1,byrow = T)
  
  # DiagSplit is a function to derive the diagram for the DaC.
  Diag_split = matrix(mapply(X_split=X_split, i1=i1, j1=j1, FUN=DiagSplit, m=m,SIMPLIFY = F),
                      m-1, m-1)
  #####################

  return(Diag_split)
}