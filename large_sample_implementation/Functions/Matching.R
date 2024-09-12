# Matching two groups of points
# Check if every points from X is in the convex hull of Y and so is Y.

Matching <- function(X,Y,eps){
  label=1
  nX=NROW(X)
  nY=NROW(Y)
  
  m=NCOL(X)
  
  #Check if X[i,] is in the convex hull of Y given eps error.
  for(i in 1:nX){
    Xi=X[i,]
    if(any(apply(Y,
                 1, function(x, want) isTRUE(all.equal(x, want)), Xi))){
      next
    }
    Yeps=Xi+(Y-Xi)*(1-eps/rowNorms(Y-Xi))
    ch=convhulln(Yeps)
    if(!inhulln(ch,matrix(Xi,1,2))){
      return(0)
    }
  }
  
  for (j in 1:nY) {
    Yj=Y[j,]
    if(any(apply(X,
             1, function(x, want) isTRUE(all.equal(x, want)), Yj))){
      next
    }
    Xeps=Yj+(X-Yj)*(1-eps/rowNorms(X-Yj))
    ch=convhulln(Xeps)
    if(!inhulln(ch,matrix(Yj,1,2))){
      return(0)
    }
  }
  
  return(label)
}