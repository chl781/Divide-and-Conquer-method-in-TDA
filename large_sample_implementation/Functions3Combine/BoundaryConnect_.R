BoundaryConnect_ <- function(i1,j1,k1,n1,dist1,dist2,dist3){
  # i1,j1,k1 is the index of the sub-region.
  # dist1, dist2, dist3 is the length of the sub-region.
  # n1 is the number of total sub-region.
  
  if(i1==1|i1==n1-1){
    if(j1==1|j1==n1-1){
      if(k1==1|k1==n1-1){
        return(matrix(c(0,0,0,
                        0,0,0,
                        0,0,0),nrow=3))
      } else{
        return(matrix(c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,dist3,
                        0,0,dist3,0),nrow=4))
      }
    } else{
      if(k1==1|k1==n1-1){
        return(matrix(c(0,0,0,0,
                        0,0,dist2,0,
                        0,dist2,0,0,
                        0,0,0,0),nrow=4))
      } else{
        return(matrix(c(0,0,0,0,0,
                        0,0,dist2,0,0,
                        0,dist2,0,0,0,
                        0,0,0,0,dist3,
                        0,0,0,dist3,0),nrow=5))
      }
    }
  }else{
    if(j1==1|j1==n1-1){
      if(k1==1|k1==n1-1){
        return(matrix(c(0,dist1,0,0,
                        dist1,0,0,0,
                        0,0,0,0,
                        0,0,0,0),nrow=4))
      } else{
        return(matrix(c(0,dist1,0,0,0,
                        dist1,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,dist3,
                        0,0,0,dist3,0),nrow=5))
      }
    } else{
      if(k1==1|k1==n1-1){
        return(matrix(c(0,dist1,0,0,0,
                        dist1,0,0,0,0,
                        0,0,0,dist2,0,
                        0,0,dist2,0,0,
                        0,0,0,0,0),nrow=5))
      } else{
        return(matrix(c(0,dist1,0,0,0,0,
                        dist1,0,0,0,0,0,
                        0,0,0,dist2,0,0,
                        0,0,dist2,0,0,0,
                        0,0,0,0,0,dist3,
                        0,0,0,0,dist3,0),nrow=6))
      }
    }
  }
}
