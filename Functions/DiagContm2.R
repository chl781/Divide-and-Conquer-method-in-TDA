

DiagContm2 <- function(X,m,maxscale,maxdimension,range){
  #Split data into (m-1)*(m-1) parts.
  gap1=seq(range[1,1],range[1,2],length.out =m)
  gap2=seq(range[2,1],range[2,2],length.out =m)
  
  #
  X_split=matrix(list(),m-1,m-1)
  for (i in 1:(m-1)) {
    for(j in 1:(m-1)){
      X_split[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
    }
  }
  
  ##### without running with Full data
  
  #DiagRips <- ripsDiag(
  #  X = X, maxdimension = maxdimension, maxscale = maxscale,
  #  library = "Dionysus", location = TRUE, printProgress = F)
  #if( any(DiagRips$diagram[,1]==1) ){
  #  index=which(DiagRips$diagram[,1]==1)
  #  t=index[which.max(DiagRips$diagram[index,3]-DiagRips$diagram[index,2])]
  #  S[1,]=c(DiagRips$diagram[t,2],DiagRips$diagram[t,3])
  #}else{
  #  S[1,]=c(0,0)
  #}
  
  ####### Begin to split the data into parts.
  
  Diag_split=matrix(list(),m-1,m-1)
   
  for(i1 in 1:(m-1)){
    for (j1 in 1:(m-1)) {
      X_split_ij=X_split[[i1,j1]]
      n1=NROW(X_split_ij)
      if(n1<=3){
        next()
      }
      
      # X_split_ij
      ## Data generation
      side_number_eliminated=0 # the outer boundary will not be counted.
      if(i1==1|i1==m-1){
        side_number_eliminated=side_number_eliminated+1
      }
      if(j1==1|j1==m-1){
        side_number_eliminated=side_number_eliminated+1
      }
      
      dist1=matrix(0,n1+4-side_number_eliminated,n1+4-side_number_eliminated)
      for (i in 1:n1) {
        for(j in i:n1){
          # Triangle matrix
          dist1[i,j]=dist1[j,i]=norm(X_split_ij[i,]-X_split_ij[j,],type="2")
        }
      }
      
      i=n1+1
      if(i1!=1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(X_split_ij[j,1]-gap1[i1])
        }
        i=i+1
      }
      
      if(i1!=m-1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(gap1[i1+1]-X_split_ij[j,1])
        }
        i=i+1
      }
      
      if(j1!=1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(X_split_ij[j,2]-gap2[j1])
        }
        i=i+1
      }
      
      if(j1!=m-1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(gap2[j1+1]-X_split_ij[j,2])
        }
      }
      
      ### Set the connection between different boundaries.
      
      dist1[(n1+1):(n1+4-side_number_eliminated),
            (n1+1):(n1+4-side_number_eliminated)]=
        BoundaryConnect(i1,j1,m,gap1[i1+1]-gap1[i1],gap2[j1+1]-gap2[j1])
      
      
      # Change X to be distance matrix
      Diag_split[[i1,j1]] <- ripsDiag(
        X = dist1, maxdimension = maxdimension, maxscale = maxscale,
        library = "Dionysus",dist="arbitrary", location = TRUE, printProgress = F)
    }
  }
  #####################

  return(Diag_split)
}