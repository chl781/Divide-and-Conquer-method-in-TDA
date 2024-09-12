# This file is served as utilities functions.

##### The following is to derive the ind_suspicious and ind_subfeature

ind_find <- function(Diag_split,X_split,m) {
  ind_suspicious_ij = c() # temporary variable for ind_suspicious
  ind_subfeature_ij = c() # temporary variable for ind_subfeature
  length_subfeature = c() # Record the length of subfeature feature
  length_suspicious = c() # Record the length of suspicious feature
  
  lengthX = NROW(X_split)
  Diag_split_ij = Diag_split
  index = which(Diag_split_ij$diagram[,1]==1)
  
  if(length(index) == 0){
    return(list(NULL,NULL,0,0))
  }
  for (ij in 1:length(index)){
    ########### Check if there is a boundary contained in the feature.
    if((lengthX+1) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+2) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+3) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])|
       (lengthX+4) %in% as.vector(Diag_split_ij[["cycleLocation"]][[ index[ij] ]])){
      ###################
      
      # Rule out the case of ij if birth time is 0. 
      # This means that this is fake feature which is composed only by the boundaries.
      if(Diag_split_ij$diagram[index[ij],2] != 0){
        ind_suspicious_ij = c(ind_suspicious_ij, index[ij])
      }
      ###################
    }else{
      ind_subfeature_ij = c(ind_subfeature_ij, index[ij])
    }
  }
  
  length_subfeature = length(ind_subfeature_ij)
  length_suspicious = length(ind_suspicious_ij)
  
  return(list(ind_suspicious_ij,ind_subfeature_ij,length_subfeature,length_suspicious))
}


##### Subfeature PD finding.

# i is used to record the current position.
subfeature_Diag <- function(Diag_split, ind_subfeature){
  birth = as.numeric(Diag_split$diagram[ind_subfeature,2])
  death = as.numeric(Diag_split$diagram[ind_subfeature,3])
  return(list(birth,death))
}


##### Retrieve the data from the PD for sub-feature data

## Utils
extract <- function(Y,X){
  return(X[as.vector(unique(unlist(Y))),])
}

unlist_part <- function(X){
  t=list()
  for(i in 1:length(X)){
    if(is.list(X[[i]])){
      t=c(t,X[[i]])
    }else{
      t=c(t,X[i])
    }
  }
  return(t)
}

# Get rid of the short list.
drop_list <- function(X){
  if(NROW(X)>=4){
    return(X)
  }
}

########################################

subfeature_cycle_Find <- function(X_split,Diag_split, ind_subfeature){
  #subfeature_cycles = list()
  if(is.null(ind_subfeature)){
    return(NULL)
  }
  subfeature_cycles = Diag_split$cycleLocation[ind_subfeature]
  X_subfeature = mapply(Y=subfeature_cycles,FUN=extract,X=rep(list(X_split),
                                            length(subfeature_cycles)),SIMPLIFY = F )
  return(X_subfeature)
}


############### 
############### The following is suspicious features

suspicious_Find <- function(X_split,Diag_split,ind_suspicious){
  if(is.null(ind_suspicious)){
    return(NULL)
  }
  l1=NROW(X_split)
  t=Diag_split$cycleLocation[ind_suspicious]
  #print(t)
  t=mapply(t,num = l1, FUN = is.small)
  if(is.list(t)){
    X_suspicious=mapply(Y=t,FUN=extract,X=rep(list(X_split),
                                              length(t)),SIMPLIFY = F )
  } else{
    X_suspicious=X_split[t,]
  }
  
  #print(X_suspicious)
  return(X_suspicious)
}

is.small <- function(t,num){
  return(t[t<=num])
}

############# Below is to construct the dist matrix.

dist_construct <- function(X,Y){
  dist1=dist(X,Y)
  return(min(dist1))
}
