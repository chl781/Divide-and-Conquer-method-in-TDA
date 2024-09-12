# Matching two groups of points
# Check if every points from X is in the convex hull of Y and so is Y.
# For this version code, we only consider if there are points close enough from the other group.


Matching1 <- function(X,Y,eps){
  label=1
  nX=NROW(X)
  nY=NROW(Y)
  
  #m=NCOL(X)
  
  #Check if X[i,] is in the convex hull of Y given eps error.
  for(i in 1:nX){
    Xi=X[i,]
    if(any(apply(Y,
                 1, function(x, want) isTRUE(all.equal(x, want)), Xi))){
      next
    }
    if(any(rowNorms(Y-matrix(rep(Xi,NROW(Y)),
                             ncol=2,byrow = T))>error)){
      return(0)
    }
  }
  
  for (j in 1:nY) {
    Yj=Y[j,]
    if(any(apply(X,
                 1, function(x, want) isTRUE(all.equal(x, want)), Yj))){
      next
    }
    if(any(rowNorms(X-matrix(rep(Yj,NROW(X)),
                             ncol=2,byrow = T))>error)){
      return(0)
    }
  }
  
  return(label)
}