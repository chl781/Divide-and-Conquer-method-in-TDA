
# We introduce the boundary and combine the features within each sub-region.
Diagm3Combine_Withboundary_Combine_Subregion_ <- function(X_split,m,Diag_split,
                                        range,maxdimension,maxscale,error){
  # Find the ind_suspicious and ind_subfeature
  
  ind1 = mapply(Diag_split = Diag_split, X_split = X_split,
                FUN = ind_find_Withboundary_, m = m, SIMPLIFY = F)
  
  ind_suspicious = array(sapply(ind1,"[[",1),c(m-1,m-1,m-1))
  ind_subfeature = array(sapply(ind1,"[[",2),c(m-1,m-1,m-1))
  length_subfeature = array(sapply(ind1,"[[",3),c(m-1,m-1,m-1)) # This is the length for the entries in the ind_subfeature
  length_bound_matrix = array(sapply(ind1,"[[",4),c(m-1,m-1,m-1)) # This is the length for the ind_suspicious.
  
  length_bound = sum(length_bound_matrix) # length_bound: the total length of ind_suspicious
  
  ##################### derive the PD for the sub-features data
  
  subfeature_PD = matrix(0, ncol = 3, nrow = sum(length_subfeature))
  colnames(subfeature_PD) <- c("dimension","birth","death")
  subfeature_PD[,1] = 2 # Only records for dimension=2
  if(sum(length_subfeature) > 0){
    bd = mapply(Diag_split = Diag_split, ind_subfeature = ind_subfeature,
                FUN = subfeature_Diag, SIMPLIFY = F)
    subfeature_PD[,2] = unlist(sapply(bd,"[[",1))
    subfeature_PD[,3] = unlist(sapply(bd,"[[",2))
  }
  
  # Retrieve the data from the PD for sub-feature data
  subfeature_cycles = vector("list",sum(length_subfeature))
  
  if(sum(length_subfeature) > 0){
    subfeature_cycles = mapply(X_split = X_split,Diag_split = Diag_split, 
                               ind_subfeature = ind_subfeature,
                               FUN = subfeature_cycle_Find, SIMPLIFY = F)
    subfeature_cycles = unlist_part(subfeature_cycles[!sapply(subfeature_cycles,is.null)])
  }
  
  
  ###################### The following is to merge the suspicious data
  
  
  ##
  
  # extract all of cycles data 
  X_suspicious=list()
  num=1
  for(i in 1:(m-1)){
    for(j in 1:(m-1)){
      for(k in 1:(m-1)){
        if(length_bound_matrix[i,j,k]>0){
          for(ij in 1:length_bound_matrix[i,j,k]){
            index11=unique(as.vector(Diag_split[[i,j,k]]$cycleLocation[[ind_suspicious[[i,j,k]][ij]]]))
            index11=index11[index11<=NROW(X_split[[i,j,k]])]
            if(length(X_suspicious)<num){
              X_suspicious[[num]]=X_split[[i,j,k]][index11,]
            }else{
              X_suspicious[[num]]=rbind(X_suspicious[[num]],X_split[[i,j,k]][index11,])
            }
          }
          num=num+1
        }
      }
    }
  }
  
  ## +6 means the extra boundary.
  length_bound=length(X_suspicious)
  if(length_bound>0){
    dist_bound1 = matrix(Inf,length_bound+6,length_bound+6) # Construct the distance matrix among all of the suspicious features.
    diag(dist_bound1)=0 # Denote the distance of feature itself is 0
  }
  
  for(i in 1:(length_bound-1)){
    for(j in (i+1):length_bound){
      dist_bound1[j,i]=dist_bound1[i,j]=dist_construct(X_suspicious[[i]],X_suspicious[[j]])
    }
  }
  
  i=1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,1]-gap1[1]))
    if(min(abs(X_suspicious[[j]][,1]-gap1[1]))==Inf){
      print(c(i,j))
    }
  }
  
  i=i+1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,1]-gap1[m]))
  }
  
  i=i+1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,2]-gap2[1]))
  }
  
  i=i+1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,2]-gap2[m]))
  }
  
  
  i=i+1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,3]-gap3[1]))
    if(min(abs(X_suspicious[[j]][,1]-gap3[1]))==Inf){
      print(c(i,j))
    }
  }
  
  i=i+1
  for(j in 1:length_bound){
    dist_bound1[length_bound+i,j]=dist_bound1[j,length_bound+i]=min(abs(X_suspicious[[j]][,3]-gap3[m]))
  }
  
  dist_bound1[(length_bound+1):(length_bound+6),(length_bound+1):(length_bound+6)]=BoundaryConnect_(2, 2, 2, m, gap1[m]-gap1[1], gap2[m]-gap2[1], gap3[m]-gap3[1])
  
  # Save the distance matrix
  # write.table(dist_bound1,file = "dist_bound.txt",row.names = F,col.names = F,sep=",")
  
  ################
  
  # Notice that the maxscale in the following function cannot be too large.
  # Running the Rips filtration based on the dist matrix among different suspicious features.
  
  diag_suspicious = ripsDiag(dist_bound, maxdimension, maxscale = maxscale,
                             library = "Dionysus", dist="arbitrary",
                             location = T)
  suspicious_ind = which(diag_suspicious$diagram[,1] == 2) # Only allow dimension to be 2
  Combined = vector("list",length(suspicious_ind))
  Combined_diag_indices = vector("list",length(suspicious_ind))
  Combined_bound = vector("list",length(suspicious_ind))
  
  # Record the indexes that save the data in suspicious data.
  non_empty = c()
  for(i1 in 1:(m-1)){ #j1 here
    for(j1 in 1:(m-1)){ #i1 here
      for(k1 in 1:(m-1)){
        t=1
        if(length_bound_matrix[i1,j1,k1]!=0){
          non_empty = rbind(non_empty,
                            cbind(matrix(rep(c(i1,j1,k1),length_bound_matrix[i1,j1,k1]),
                                         nrow = length_bound_matrix[i1,j1,k1],
                                         byrow = T),1:length_bound_matrix[i1,j1,k1]))
        }
      }
    }
  }
  
  ##### Retrieve the data
  
  
  ################# 
  #The below is trying to do for loop to retrieve the data splited for distance method.
  
  num=1
  for(one in suspicious_ind){
    Combined_ind = unique(as.vector(diag_suspicious$cycleLocation[[one]]))
    #col = 1
    for(i in Combined_ind){
      df = Diag_split[[non_empty[i,1],non_empty[i,2],non_empty[i,3]]][["cycleLocation"]] # Retrieve the index
      df1 = df[[ind_suspicious[[non_empty[i,1], non_empty[i,2], non_empty[i,3] ]][ non_empty[i,4] ] ]] 
      df1 = unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i,1], non_empty[i,2],
                                                     non_empty[i,3]]])])
      if(length(df1) <= 1){
        stop()
      }
      Suspicious_i = cbind(X_split[[non_empty[i,1], non_empty[i,2], non_empty[i,3]]][df1,]
                           #,rep(col,NROW(X_split[[non_empty[i,1],non_empty[i,2]]][df1,]))
      ) # Retrieve the data
      
      #colnames(Suspicious_i) <- c("x","y","col")
      colnames(Suspicious_i) <- c("x","y","z")
      Combined[[num]] = rbind(Combined[[num]],Suspicious_i)
      Combined_diag_indices[[num]] = rbind(Combined_diag_indices[[num]],
                                           non_empty[i,]) # This saves all of points constructing the loops.
      #col=col+1
    }
    num = num+1
  }
  
  ##############  
  
  ##### Comment out the projected method
  
  
  # Combine the lists
  Combined1 <- Combined
  Combined_diag_indices1=Combined_diag_indices
  t<-c()
  if(length(Combined1)>0){
    for(i in 1:length(Combined1)){
      if(is.null(Combined1[[i]])){
        t = c(t,i)
      }
    }
  }
  if(length(t)>0){
    Combined1=Combined1[-t]
    Combined_diag_indices1=Combined_diag_indices1[-t]
  }
  
  #########
  # Recover the suspicious estimates below.
  
  
  PD2 <- matrix(0,ncol=3,nrow=length(Combined1))
  colnames(PD2) <- c("dimension","birth","death")
  PD2[,1]=2 # Only allow dimension to be 2
  for (i in 1:NROW(PD2) ) {
    PD2[i,2] = BirthRecal2_(Diag_split,Combined_diag_indices[[i]],ind_suspicious)
    PD2[i,3] = DeathRecal0_(unique(Combined1[[i]])[,c(1,2,3)]) # This is a specific death estimate for d=2.
    print(i)
    # DeathRecal0 uses 4 for loops, which can definitely be improved.
  }
  PD2
  
  # Combine the combined PD and PD for sub-features.
  PD_=rbind(subfeature_PD,PD)
  # The following is the representative points.
  cycle_points=append(subfeature_cycles,Combined1)
  
  ################ Second combination step.
  
  # Extract all of representative points
  RData=vector("list",length = 37)
  
  num=1
  for(i1 in 1:(m-1)){ #j1 here
    for(j1 in 1:(m-1)){ #i1 here
      for(k1 in 1:(m-1)){
        
      }
    }
  }
  # Construct the distance matrix and Run the rips filtration
  
  # Derive the PD
  
  return(list(diagram=PD_,cycle_points))
}