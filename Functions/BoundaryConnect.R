# This file is to build the connection status among the boundaries.

BoundaryConnect <- function(i1,j1,n1,dist1,dist2){
  if(i1==1|i1==n1-1){
    if(j1==1|j1==n1-1){
      return(matrix(c(0,0,0,0),nrow=2))
    } else{
      return(matrix(c(0,0,0,0,0,dist2,0,dist2,0),nrow=3))
    }
  }else{
    if(j1==1|j1==n1-1){
      return(matrix(c(0,dist1,0,dist1,0,0,0,0,0),nrow=3))
    } else{
      return(matrix(c(0,dist1,0,0,dist1,0,0,0,0,0,0,dist2,0,0,dist2,0),nrow=4))
    }
  }
}
