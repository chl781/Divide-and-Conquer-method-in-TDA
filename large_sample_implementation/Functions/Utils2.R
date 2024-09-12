# This is the utils file for DiagContm2 for parallel purpose
# po records the position of current X_split

DiagSplit <- function(X_split_ij,i1,j1,m){
  #X_split_ij = X_split[[i1,j1]]
  length_X_split_ij = NROW(X_split_ij)
  
  # Prevent trivial case
  if(length_X_split_ij < 4){ 
    return(NULL)
  }
  
  ## Data generation
  side_number_eliminated = 0 # Decide how many boundaries will be taken into account.
  if(i1 == 1 | i1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  if(j1 == 1 | j1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  
  dist1 = matrix(0, length_X_split_ij+4-side_number_eliminated,
                 length_X_split_ij+4-side_number_eliminated)
  for (i in 1:(length_X_split_ij-1)) {
    for(j in (i+1):length_X_split_ij){
      # Construct distance among points.
      dist1[i,j] = dist1[j,i] = norm(X_split_ij[i,] - X_split_ij[j,], type="2")
    }
  }
  
  # i is the index for the boundary.
  i = length_X_split_ij + 1
  if(i1 != 1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(X_split_ij[j,1]-gap1[i1])
    }
    i = i+1
  }
  
  if(i1 != m-1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(gap1[i1+1]-X_split_ij[j,1])
    }
    i = i+1
  }
  
  if(j1 != 1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(X_split_ij[j,2] - gap2[j1])
    }
    i = i+1
  }
  
  if(j1 != m-1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(gap2[j1+1]-X_split_ij[j,2])
    }
  }
  
  ### Compute the distance matrix among boundaries.
  
  dist1[(length_X_split_ij+1):(length_X_split_ij + 4 - side_number_eliminated),
        (length_X_split_ij+1):(length_X_split_ij + 4 - side_number_eliminated)] = 
    BoundaryConnect(i1, j1, m, gap1[i1+1]-gap1[i1], gap2[j1+1]-gap2[j1])
  
  
  # Use distance matrix to compute the diagram.
  Diag_split_ij <- ripsDiag(
    X = dist1, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", dist="arbitrary", location = TRUE, printProgress = F)
  return(Diag_split_ij)
}

