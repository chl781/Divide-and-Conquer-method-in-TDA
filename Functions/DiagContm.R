# This code is to design for cutting into n1*n1 boxes and then combine it.

# m is the step size, choice is about how to identify the indices to combine.
DiagContm <- function(X,m,maxscale,maxdimension,range,choice){
  #Split data into 4 parts.
  A1=range[1,1]
  A2=range[1,2]
  B1=range[2,1]
  B2=range[2,2]
  gap1=seq(A1,A2,length.out =m)
  gap2=seq(B1,B2,length.out =m)
  
  X1=matrix(list(),m-1,m-1)
  for (i in 1:(m-1)) {
    for(j in 1:(m-1)){
      X1[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
    }
  }
  
  S=matrix(0,nrow = 4,ncol = 2)
  
  ##### Full data
  
  DiagRips <- ripsDiag(
    X = X, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips$diagram[,1]==1) ){
    index=which(DiagRips$diagram[,1]==1)
    t=index[which.max(DiagRips$diagram[index,3]-DiagRips$diagram[index,2])]
    S[1,]=c(DiagRips$diagram[t,2],DiagRips$diagram[t,3])
  }else{
    S[1,]=c(0,0)
  }
  
  ####### Begin to split the data into parts.
  
  Diag2=matrix(list(),m-1,m-1)
  t2=matrix(0,m-1,m-1)
  
  for(i1 in 1:(m-1)){
    for (j1 in 1:(m-1)) {
      X2=X1[[i1,j1]]
      n1=NROW(X2)
      if(n1<=2){
        next()
      }
      
      # X2
      ## Data generation
      si=sj=0 # the very boundary will not be counted.
      if(i1==1|i1==m-1){
        si=1
      }
      if(j1==1|j1==m-1){
        sj=1
      }
      
      dist1=matrix(0,n1+4-si-sj,n1+4-si-sj)
      for (i in 1:n1) {
        for(j in i:n1){
          dist1[i,j]=dist1[j,i]=norm(X2[i,]-X2[j,],type="2")
        }
      }
      
      i=n1+1
      if(i1!=1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(X2[j,1]-gap1[i1])
        }
        i=i+1
      }
      
      if(i1!=m-1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(gap1[i1+1]-X2[j,1])
        }
        i=i+1
      }
      
      if(j1!=1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(X2[j,2]-gap2[j1])
        }
        i=i+1
      }
      
      if(j1!=m-1){
        for(j in 1:n1){
          dist1[i,j]=dist1[j,i]=abs(gap2[j1+1]-X2[j,2])
        }
      }
      
      ### Set the connection between different boundaries.
      
      dist1[(n1+1):(n1+4-si-sj),(n1+1):(n1+4-si-sj)]=BoundaryConnect(i1,j1,m)
      
      
      # Change X to be distance matrix
      Diag2[[i1,j1]] <- ripsDiag(
        X = dist1, maxdimension = maxdimension, maxscale = maxscale,
        library = "Dionysus",dist="arbitrary", location = TRUE, printProgress = F)
      # Be careful of the combination step
      DiagRips2=Diag2[[i1,j1]]
      if( any(DiagRips2$diagram[,1]==1) ){
        index=which(DiagRips2$diagram[,1]==1)# Change it to dimension = 1.
        t2[i1,j1]=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
      }else{
        t2[i1,j1]=-1
      }
    }
  }
  #####################
  # Combine the data. 
  # Notice that this step still not involved in the combination step.
  # 1st method
  if(choice==1){
    ind1=c()
    for (i1 in 1:(m-1)) {
      for (j1 in 1:(m-1)) {
        if(t2[i1,j1]>0){
          ind1=rbind(ind1,c(i1,j1))
        }
      }
    }
    
    represent=ripsDiag(
      X = ind1, maxdimension = maxdimension, maxscale = m,
      library = "Dionysus", location = TRUE, printProgress = F)
    
    if( any(represent$diagram[,1]==1) ){
      index=which(represent$diagram[,1]==1)# Change it to dimension = 1.
      trepresent=index[which.max(represent$diagram[index,3]-represent$diagram[index,2])]
    }else{
      trepresent=-1
    }
    
    indexComb1=matrix(as.vector(represent[["cycleLocation"]][[trepresent]]),ncol=2,byrow = F)
  }
  # indexComb1 is the index sets that contain the features.
  
  
  ### 2nd method
  # Combining method based on the boundary 
  # Boundary should be saved in the previous procedure.
  if(choice==2){
    # Extract the boundary information.
    bound1=BoundaryDetect(t2,X1,Diag2,m)
  }
  
  ########
  
  # Run the combining programming.
  if(choice==1){
    data=c()
    for(i1 in 1:(m-1)){
      for(j1 in 1:(m-1)){
        if(t2[i1,j1]!=-1){
          X2=rbind(X1[[i1,j1]],c(0,0),c(0,0),c(0,0),c(0,0))
          data1=unique(X2[as.vector(Diag2[[i1,j1]][["cycleLocation"]][[t2[i1,j1]]]),])
          data=rbind(data,data1)
        }
      }
    }
  } 
  
  # Combining the data based on the second method
  if(choice==2){
    # Denote a choosing set
    ind1=matrix(1,m-1,m-1)
    ind1(t2<=0)=0
    for(i1 in 1:(m-1)){
      for(j1 in 1:(m-1)){
        if(ind1[i1,j1]==1){
          
        }
      }
    }
  }
  
  # Blindly combine. Not realistic.
  if(choice==3){
    data=c()
    for(i1 in 1:(m-1)){
      for(j1 in 1:(m-1)){
        if(j1 %in% indexComb1[indexComb1[,1]==i1,2]){
          X2=rbind(X1[[i1,j1]],c(0,0),c(0,0),c(0,0),c(0,0))
          data1=unique(X2[as.vector(Diag2[[i1,j1]][["cycleLocation"]][[t2[i1,j1]]]),])
          data=rbind(data,data1)
        }
      }
    }
  }
  
 
  
  # Here is a problem, because if there is no complete half of the data, then we cannot
  # use the sophisticated method to produce it.
  
  # Notice that this procedure do not check for the coincidence for the representative points.
  data=data[data[,1]!=0 & data[,2]!=0,]
  
  # Based on the new data, we analyze the estimated birth time and death time.
  # Re-estimate.
  DiagRips5 <- ripsDiag(
    X = data, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips5$diagram[,1]==1) ){
    index=which(DiagRips5$diagram[,1]==1)
    t5=index[which.max(DiagRips5$diagram[index,3]-DiagRips5$diagram[index,2])]
    S[2,]=c(DiagRips5$diagram[t5,2],DiagRips5$diagram[t5,3])
  }else{
    S[2,]=c(0,0)
  }
  
  # Added estimates.
  
  
  # Birth estimate is unrealistic because there could be missing sub-boxes.
  #  Too conservative for distribution indpendent case, but can fix the distribution case.
  # But, if there is sparse case, then this will decrease the time estimate.
  
  BirthRecal=0 # Initialization for BirthRecal.
  
  if(choice==1){
    for(i1 in 1:(m-1)){
      for(j1 in 1:(m-1)){
        if(t2[i1,j1]==0|t2[i1,j1]==-1){
          next()
        }
        if(BirthRecal<Diag2[[i1,j1]]$diagram[t2[i1,j1],2]){
          BirthRecal=Diag2[[i1,j1]]$diagram[t2[i1,j1],2]
        }
      }
    }
  }
  
  if(choice==2){
    for(i1 in 1:(m-1)){
      for(j1 in 1:(m-1)){
        if(j1 %in% indexComb1[indexComb1[,1]==i1,2]){
          if(BirthRecal<Diag2[[i1,j1]]$diagram[t2[i1,j1],2]){
            BirthRecal=Diag2[[i1,j1]]$diagram[t2[i1,j1],2]
          }
        }
      }
    }
  }
  
  # This birth estimate is for general.
  # Should do it in a more sophisticated way.
  Deathrecal=DeathRecal_Circle(data,1,2,3) 
  # Another approximation is to sample several points from different sub-boxes.
  S[3,]=c(BirthRecal,Deathrecal)
  
  # Approximate estimate based on whole data
  gmra = gmra.create.ipca(X, eps=0, dim=1, maxKids=1, stop=4)
  res <- multiscale.rips(gmra, maxD = 1)
  
  t3<-which.max(res$diagram[res$diagram[,4]==1,2]- 
                  res$diagram[res$diagram[,4]==1,1])
  if(sum(res$diagram[res$diagram[,4]==1,])==0){
    MSE_Linear_birth=0
    MSE_Linear_death=0
  }else{
    MSE_Linear_birth=res$diagram[which(res$diagram[,4]==1)[t3],1]
    MSE_Linear_death=res$diagram[which(res$diagram[,4]==1)[t3],2]
  }
  
  S[4,]=c(MSE_Linear_birth,MSE_Linear_death)
  
  # Check if we can combine the sub-features.
  return(S)
}