# This function serves as Combine method v3.

# Inputs:
# X_split is the data separate into several subregions.
# m is the number of subregions to divide.
# Diag_split is the persistent diagram with cycles into different subregions.
# range is the range of data that we consider.
# maxdimension and maxscale are inputs to ripsDiag
# error is te error bound allowed to use in the cancellation method.

Diagm3Combine <- function(X_split,m,Diag_split,
                          range,maxdimension,maxscale,error){
  gap1 = seq(range[1,1], range[1,2], length.out = m)
  gap2 = seq(range[2,1], range[2,2], length.out = m)
  
  # Find the ind_suspicious and ind_subfeature
  
  ind1 = mapply(Diag_split = Diag_split, X_split = X_split,
                FUN = ind_find, m = m, SIMPLIFY = F)
  
  ind_suspicious = matrix(sapply(ind1,"[[",1),m-1,m-1)
  ind_subfeature = matrix(sapply(ind1,"[[",2),m-1,m-1)
  length_subfeature = matrix(sapply(ind1,"[[",3),m-1,m-1) # This is the length for the entries in the ind_subfeature
  length_bound_matrix = matrix(sapply(ind1,"[[",4),m-1,m-1) # This is the length for the ind_suspicious.
  
  length_bound = sum(length_bound_matrix) # length_bound: the total length of ind_suspicious
  
  
  
  # for (i in 1:(m-1)) {
  #   for (j in 1:(m-1)) {
  #     ind_suspicious_ij = c() # temporary variable for ind_suspicious
  #     ind_subfeature_ij = c() # temporary variable for ind_subfeature
  #     lengthX = NROW(X_split[[i,j]])
  #     Diag_split_ij = Diag_split[[i,j]]
  #     index = which(Diag_split_ij$diagram[,1]==1)
  #     if(length(index) == 0){
  #       next
  #     }
  #     for (ij in 1:length(index)){
  #       ########### Check if there is a boundary contained in the feature.
  #       if((lengthX+1) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
  #          (lengthX+2) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
  #          (lengthX+3) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
  #          (lengthX+4) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])){
  #         ###################
  #         
  #         # Rule out the case of ij if birth time is 0. 
  #         # This means that this is fake feature which is composed only by the boundaries.
  #         if(Diag_split_ij$diagram[index[ij],2] != 0){
  #           ind_suspicious_ij = c(ind_suspicious_ij, index[ij])
  #         }
  #         ###################
  #       }else{
  #         ind_subfeature_ij = c(ind_subfeature_ij, index[ij])
  #       }
  #     }
  #     if(!is.null(ind_suspicious_ij)){
  #       ind_suspicious[[i,j]] = ind_suspicious_ij
  #     }
  #     if(!is.null(ind_subfeature_ij)){
  #       ind_subfeature[[i,j]] = ind_subfeature_ij
  #     }
  #   }
  # }
  
  
  
  ##################### derive the PD for the sub-features data
  
  subfeature_PD = matrix(0, ncol = 3, nrow = sum(length_subfeature))
  colnames(subfeature_PD) <- c("dimension","birth","death")
  subfeature_PD[,1]=1
  if(sum(length_subfeature) > 0){
    bd = mapply(Diag_split = Diag_split, ind_subfeature = ind_subfeature,
                FUN = subfeature_Diag, SIMPLIFY = F)
    subfeature_PD[,2] = unlist(sapply(bd,"[[",1))
    subfeature_PD[,3] = unlist(sapply(bd,"[[",2))
  }
  
  # Retrieve the data from the PD for sub-feature data
  subfeature_cycles = vector("list",sum(length_subfeature))
  
  if(sum(length_subfeature) > 0){
    #for(i in 1:(m-1)^2){
        # for(k in sequence(length_subfeature[i])){
        #   if(i==1){
        #     subfeature_cycles[[k]] = 
        #       Diag_split[[i]]$cycleLocation[[ ind_subfeature[[i]][k] ]]
        #   }else{
        #     subfeature_cycles[[ sum(length_subfeature[1:(i-1)])+k ]]=
        #       Diag_split[[i]]$cycleLocation[[ ind_subfeature[[i]][k] ]]
        #   }
        # }
      subfeature_cycles = mapply(X_split = X_split,Diag_split = Diag_split, 
                                 ind_subfeature = ind_subfeature,
                                 FUN = subfeature_cycle_Find, SIMPLIFY = F)
      subfeature_cycles = unlist_part(subfeature_cycles[!sapply(subfeature_cycles,is.null)])
    #}
  }
  
  
  ###################### The following is to merge the suspicious data
  
  ######### Combine the data from sub-regions by using the suspicious data but not projected boundary.
  if(sum(ind_suspicious)>0){
    # put everything there.
    dist_bound = matrix(Inf,length_bound,length_bound) # Construct the distance matrix among all of the suspicious features.
    diag(dist_bound)=0 # Denote the distance of feature itself is 0
  }
  
  ##
  
  # extract all of cycles data 
  X_suspicious <- mapply(Diag_split = Diag_split, ind_suspicious = ind_suspicious, X_split = X_split,
                         FUN = suspicious_Find, SIMPLIFY = F)
  X_suspicious = unlist_part(X_suspicious[!sapply(X_suspicious,is.null)])
  X_suspicious = mapply(X_suspicious, FUN = drop_list,SIMPLIFY = F) # drop short list
  X_suspicious = unlist_part(X_suspicious[!sapply(X_suspicious,is.null)])
  
  ##
  ## dist construction
  length_dist_bound = length(X_suspicious)
  dist_bound = matrix(mapply(FUN = dist_construct, X = rep(X_suspicious,times=length_dist_bound),
                           Y = rep(X_suspicious,each=length_dist_bound)),
                    length_dist_bound,length_dist_bound)
  
################
  ### Below is the classical for loop way to construct the dist matrix. Ignore
  
  # for (i1 in 1:((m-1)^2-1) ) { 
  #   for (j1 in (i1+1):((m-1)^2) ) { # i1, j1 are the index for bound_matrix
  #     #for(j1 in c(i1+1,i1+m-1)) {
  #     
  #     if ( length_bound_matrix[i1] == 0 | length_bound_matrix[j1] == 0 ){
  #       next
  #     }
  #     
  #     for (i2 in 1:length_bound_matrix[i1]) {
  #       for (j2 in 1:length_bound_matrix[j1]) {
  #         Diag_split_i1_i2_location = Diag_split[[i1]]$cycleLocation[[
  #           ind_suspicious[[i1]][i2] ]]
  #         Diag_split_i1_i2_location = 
  #           unique(Diag_split_i1_i2_location[ as.vector(Diag_split_i1_i2_location) <= 
  #                                             NROW(X_split[[i1]]) ])
  #         Diag_split_i1_i2 = X_split[[i1]][Diag_split_i1_i2_location,]
  #         if( NROW(Diag_split_i1_i2)<=1 | NROW(Diag_split_i1_i2) != 
  #             length(Diag_split_i1_i2_location) ){
  #           #stop()
  #           warning()
  #         }
  #         Diag_split_j1_location = Diag_split[[j1]]$cycleLocation[[ 
  #           ind_suspicious[[j1]][j2] ]]
  #         Diag_split_j1_location = 
  #           unique(Diag_split_j1_location[as.vector(Diag_split_j1_location) 
  #                                         <= NROW(X_split[[j1]])])
  #         Diag_split_j1_j2 = X_split[[j1]][Diag_split_j1_location,]
  #         if(NROW(Diag_split_j1_j2) <= 1 | NROW(Diag_split_j1_j2) != 
  #            length(Diag_split_j1_location)){
  #           #stop()
  #           warning()
  #         }
  #         
  #         # determine the position of the dist matrix.
  #         if(i1 > 1){
  #           i3 = sum(length_bound_matrix[1:(i1-1)]) + i2
  #         }
  #         if(i1 == 1){
  #           i3 = i2
  #         }
  #         j3 = sum(length_bound_matrix[1:(j1-1)]) + j2
  # 
  #         dist_bound[i3,j3] = dist_bound[j3,i3] = 
  #           min( dist(Diag_split_i1_i2,Diag_split_j1_j2) )
  #       }
  #     }
  #   }
  # }
####################
  
  ## Constuct dist_matrix for the suspicious features within each block
  # for (i1 in 1:(m-1)^2 ) {
  #   length_bound_i1 = length_bound_matrix[i1]
  #   if(length_bound_i1>1){
  #     for(i2 in 1:(length_bound_i1-1) ){
  #       for(j2 in (i2+1):length_bound_i1 ){
  #         Diag_split_i1_i2_location = Diag_split[[i1]][["cycleLocation"]][[ ind_suspicious[[i1]][i2] ]]
  #         Diag_split_i1_i2_location = unique(Diag_split_i1_i2_location[as.vector(Diag_split_i1_i2_location)<=NROW(X_split[[i1]])])
  #         Diag_split_i1_i2 = X_split[[i1]][Diag_split_i1_i2_location,] 
  #         # Diag_split_i1_i2 is the data represent for ind_suspicious[[i1]][i2]
  #         if(NROW(Diag_split_i1_i2) <= 1 | 
  #            NROW(Diag_split_i1_i2) != length(Diag_split_i1_i2_location)){
  #           stop()
  #         }
  #         
  #         Diag_split_i1_j2_location = 
  #           Diag_split[[i1]][["cycleLocation"]][[ ind_suspicious[[i1]][j2] ]]
  #         Diag_split_i1_j2_location = unique(Diag_split_i1_j2_location[
  #           as.vector(Diag_split_i1_j2_location)<=NROW(X_split[[i1]])])
  #         df22 = X_split[[i1]][Diag_split_i1_j2_location,] # df22 is the data represent for ind_suspicious[[i1]][j2]
  #         if( NROW(df22)<=1 | NROW(df22)!=length(df2) ){
  #           stop()
  #         }
  #         
  #         if(i1>1){
  #           i3 = sum(length_bound_matrix[1:(i1-1)])+i2
  #           j3 = sum(length_bound_matrix[1:(i1-1)])+j2
  #         }
  #         if(i1==1){
  #           i3 = i2
  #           j3 = j2
  #         }
  #         dist_bound[i3,j3] = min(dist(Diag_split_i1_i2,df22))
  #         dist_bound[j3,i3] = dist_bound[i3,j3]
  #       }
  #     }
  #   }
  # }
  
  # The above is the for loop way to construct the dist matrix.
################
    
  # Notice that the maxscale in the following function cannot be too large.
  # Running the Rips filtration based on the dist matrix among different suspicious features.
  
  diag_suspicious = ripsDiag(dist_bound, maxdimension, maxscale = maxscale,
                           library = "Dionysus", dist="arbitrary",
                           location = T)
  suspicious_ind = which(diag_suspicious$diagram[,1] == 1)
  Combined = vector("list",length(suspicious_ind))
  Combined_diag_indices = vector("list",length(suspicious_ind))
  Combined_bound = vector("list",length(suspicious_ind))
  
  # Record the indexes that save the data in suspicious data.
  non_empty = c()
  for(j1 in 1:(m-1)){ #j1 here
    for(i1 in 1:(m-1)){ #i1 here
      t=1
      if(length_bound_matrix[i1,j1]!=0){
        non_empty = rbind(non_empty,
                        cbind(matrix(rep(c(i1,j1),length_bound_matrix[i1,j1]),
                                     nrow = length_bound_matrix[i1,j1],
                                     byrow = T),1:length_bound_matrix[i1,j1]))
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
      df = Diag_split[[non_empty[i,1],non_empty[i,2]]][["cycleLocation"]] # Retrieve the index
      df1 = df[[ind_suspicious[[non_empty[i,1], non_empty[i,2] ]][ non_empty[i,3] ] ]] 
      df1 = unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i,1],non_empty[i,2]]])])
      if(length(df1) <= 1){
        stop()
      }
      Suspicious_i = cbind(X_split[[non_empty[i,1],non_empty[i,2]]][df1,]
                           #,rep(col,NROW(X_split[[non_empty[i,1],non_empty[i,2]]][df1,]))
                           ) # Retrieve the data
      
      #colnames(Suspicious_i) <- c("x","y","col")
      colnames(Suspicious_i) <- c("x","y")
      Combined[[num]] = rbind(Combined[[num]],Suspicious_i)
      Combined_diag_indices[[num]] = rbind(Combined_diag_indices[[num]],
                                         non_empty[i,]) # This saves all of points constructing the loops.
      #col=col+1
    }
    num = num+1
  }

  ##############  
    
  ##### The following code tries to use cancellation method to merge. Only allow using
  ##### 2 or 3 suspicious features to merge.
  
  # Create bound for suspicious features.
  Suspicious_bound = boundFind(ind_suspicious,Diag_split,X_split,m,gap1,gap2,error)
  bound=Suspicious_bound[[1]]
  Suspicious=Suspicious_bound[[2]]
  
  # The following uses the projected method.
  Projected_Merge_ <- Projected_Merge(ind_suspicious,bound,range,Diag_split,X_split,Suspicious,non_empty
                                      ,m,gap1,gap2)
  Projected_Merge1=Projected_Merge_[[1]]
  Combined_diag_indices2=Projected_Merge_[[2]]
  
  # Combine the lists
  Combined1 <- append(Combined,Projected_Merge1)
  Combined_diag_indices1=append(Combined_diag_indices,Combined_diag_indices2)
  t<-c()
  for(i in 1:length(Combined1)){
    if(is.null(Combined1[[i]])){
      t = c(t,i)
    }
  }
  if(length(t)>0){
    Combined1=Combined1[-t]
    Combined_diag_indices1=Combined_diag_indices1[-t]
  }
  
  #########
  # Recover the suspicious estimates below.
  
  
  PD <- matrix(0,ncol=3,nrow=length(Combined1))
  colnames(PD) <- c("dimension","birth","death")
  PD[,1]=1
  for (i in 1:NROW(PD) ) {
    PD[i,2] = BirthRecal2(Diag_split,Combined_diag_indices[[i]],ind_suspicious)
    PD[i,3] = DeathRecal0(unique(Combined1[[i]])[,c(1,2)]) # This is a specific death estimate.
  }
  
  # Combine the combined PD and PD for sub-features.
  PD_=rbind(subfeature_PD,PD)
  # The following is the representative points.
  cycle_points=append(subfeature_cycles,Combined1)
  
  return(list(diagram=PD_,cycle_points))
}