# Matching several groups of points
# Check if every points from one group can be expressed by the points from other groups.
# For this version code, we only consider if there are points close enough from the other group.


Matching2 <- function(list_X,eps){
  label=1
  n1=length(list_X)
  
  #Check if X[i,] is in the convex hull of Y given eps error.
  for(i in 1:n1){
    indicator=rep(1,NROW(list_X[[i]]))
    for(j in 1:NROW(list_X[[i]])){
      Xi=list_X[[i]][j,]
      for(k in 1:n1 ){
        if(k == i){
          next
        }
        if(any(is.na(list_X[[k]]-Xi)) | !any(rowNorms(list_X[[k]]-matrix(rep(Xi,NROW(list_X[[k]])),
                                                                      ncol=2,byrow = T))<=eps)){
          indicator[j]=0
        }
      }
      if(sum(indicator)<=NROW(list_X[[i]])-1){
        return(0)
      }
    }
  }
  
  return(label)
}