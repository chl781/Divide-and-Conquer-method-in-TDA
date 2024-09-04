# Distinguish the 2-d features in Diag_split as embryonic features or sub-features.

# Diag_split is Diag_split[[i1,j1,k1]] as one divided persistent diagram.
ind_find_ <- function(Diag_split,X_split,m) {
  ind_suspicious_ij = c() # temporary variable for ind_suspicious
  ind_subfeature_ij = c() # temporary variable for ind_subfeature
  length_subfeature = c() # Record the length of subfeature feature
  length_suspicious = c() # Record the length of suspicious feature
  
  lengthX = NROW(X_split)
  Diag_split_ij = Diag_split
  index = which(Diag_split_ij$diagram[,1] == 2) # Only handle the dimension = 2 case.
  
  if(length(index) == 0){
    return(list(NULL,NULL,0,0))
  }
  
  # The criterion is whether there is a supplemental boundary contained in the representative data.
  for (ij in 1:length(index)){
    ########### Check if there is a boundary contained in the feature.
    if((lengthX+1) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+2) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+3) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+4) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+5) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+6) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])){
      
      # Rule out the case of ij if birth time is 0. 
      # This means that this is fake feature which is composed only by the boundaries.
      if(Diag_split_ij$diagram[index[ij],2] != 0){
        ind_suspicious_ij = c(ind_suspicious_ij, index[ij])
      }
    }else{
      ind_subfeature_ij = c(ind_subfeature_ij, index[ij])
    }
  }
  
  length_subfeature = length(ind_subfeature_ij)
  length_suspicious = length(ind_suspicious_ij)
  
  return(list(ind_suspicious_ij,ind_subfeature_ij,length_subfeature,length_suspicious))
}


# This function changes the list within list to only list
unlist_part_ <- function(X){
  t=list()
  for(i in 1:length(X)){
    if(is.list(X[[i]])){
      t=append(t,X[[i]])
    }else{
      t=append(t,list(X[[i]]))
    }
  }
  return(t)
}