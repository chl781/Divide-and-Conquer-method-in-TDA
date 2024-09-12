# This file serves as a 2 dimensional object recovery.

DiagContm2_ <- function(X,m,maxscale,maxdimension,range){
  # Split data into (m-1)*(m-1) parts.
  gap1=seq(range[1,1], range[1,2], length.out = m)
  gap2=seq(range[2,1], range[2,2], length.out = m)
  gap3=seq(range[3,1], range[3,2], length.out = m)
  
  #
  X_split = array(list(),c(m-1,m-1,m-1))
  for (i in 1:(m-1)) {
    for(j in 1:(m-1)){
      for(k in 1:(m-1)){
        X_split[[i,j,k]] = X[ X[,1] >= gap1[i] & X[,1] < gap1[i+1] 
                            & X[,2] >= gap2[j] & X[,2] < gap2[j+1]
                            & X[,3] >= gap3[k] & X[,3] < gap3[k+1],]
      }
    }
  }
  
  ####### Begin to split the data into parts.
  
  i1 = array(rep(1:(m-1),(m-1)^2),c(m-1,m-1,m-1))
  j1 = array(rep(rep(1:(m-1),(m-1)),each=m-1),c(m-1,m-1,m-1))
  k1 = array(rep(1:(m-1),each=(m-1)^2),c(m-1,m-1,m-1))
  
  # DiagSplit is a function to derive the diagram for the DaC.
  Diag_split = array(mapply(X_split=X_split, i1=i1, j1=j1, k1=k1, FUN=DiagSplit_, m=m,SIMPLIFY = F),
                      c(m-1, m-1, m-1))
  #####################
  
  return(Diag_split)
}