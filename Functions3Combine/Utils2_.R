
# This function generate the persistent diagram for individual sub-region with supplemental boundary.

DiagSplit_ <- function(X_split_ij,i1,j1,k1,m){
  #X_split_ij = X_split[[i1,j1,k1]] is the data in sub-region. 
  length_X_split_ij = NROW(X_split_ij)
  
  # Prevent trivial case
  if(length_X_split_ij <= 4){ 
    return(NULL)
  }
  
  ## Data generation for the supplemental boundary
  side_number_eliminated = 0 # Decide how many boundaries will be taken into account.
  if(i1 == 1 | i1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  if(j1 == 1 | j1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  if(k1 == 1 | k1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  
  # Construct dist matrix among data points and supplemental boundaries. The size is determined by which sub-region is considered.
  dist1 = matrix(0, length_X_split_ij+6-side_number_eliminated,
                 length_X_split_ij+6-side_number_eliminated)
  
  # Assign weights
  for (i in 1:length_X_split_ij) {
    for(j in i:length_X_split_ij){
        # Construct distance among points.
        dist1[i,j] = dist1[j,i] = norm(X_split_ij[i,] - X_split_ij[j,], type="2")
    }
  }
  
  # Decide what is the distance between the supplementary boundary and data points.
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
    i=i+1
  }
  
  if(k1 != 1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(X_split_ij[j,3] - gap3[k1])
    }
    i = i+1
  }
  
  if(k1 != m-1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(gap3[k1+1]-X_split_ij[j,3])
    }
    i=i+1
  }
  
  ### Compute the distance matrix among boundaries.
  dist1[(length_X_split_ij+1):(length_X_split_ij + 6 - side_number_eliminated),
        (length_X_split_ij+1):(length_X_split_ij + 6 - side_number_eliminated)] = 
    BoundaryConnect_(i1, j1, k1, m, gap1[i1+1]-gap1[i1], gap2[j1+1]-gap2[j1], gap3[k1+1]-gap3[k1])
  
  
  # Use distance matrix to compute the diagram.
  Diag_split_ij <- ripsDiag(
    X = dist1, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", dist="arbitrary", location = TRUE, printProgress = F)
  return(Diag_split_ij)
}

# This function is similar to function DiagSplit_, but runs the program in a parallel fashion.
DiagSplit_parallel_ <- function(X_split_ij,i1,j1,k1,m,gap1,gap2,gap3){
  #X_split_ij = X_split[[i1,j1,k1]]
  length_X_split_ij = NROW(X_split_ij)
  
  # Prevent trivial case
  if(length_X_split_ij <= 4){ 
    return(NULL)
  }
  
  ## Decide how many supplemental boundaries there are.
  side_number_eliminated = 0 # Decide how many boundaries will be taken into account.
  if(i1 == 1 | i1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  if(j1 == 1 | j1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  if(k1 == 1 | k1 == m-1){
    side_number_eliminated = side_number_eliminated + 1
  }
  
  ## Distance matrix generation
  dist1 = matrix(0, length_X_split_ij+6-side_number_eliminated,
                 length_X_split_ij+6-side_number_eliminated)
  for (i in 1:length_X_split_ij) {
    for(j in i:length_X_split_ij){
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
    i=i+1
  }
  
  if(k1 != 1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(X_split_ij[j,3] - gap3[k1])
    }
    i = i+1
  }
  
  if(k1 != m-1){
    for(j in 1:length_X_split_ij){
      dist1[i,j] = dist1[j,i] = abs(gap3[k1+1]-X_split_ij[j,3])
    }
    i=i+1
  }
  
  ### Compute the distance matrix among boundaries.
  
  dist1[(length_X_split_ij+1):(length_X_split_ij + 6 - side_number_eliminated),
        (length_X_split_ij+1):(length_X_split_ij + 6 - side_number_eliminated)] = 
    BoundaryConnect_(i1, j1, k1, m, gap1[i1+1]-gap1[i1], gap2[j1+1]-gap2[j1], gap3[k1+1]-gap3[k1])
  
  
  # Use distance matrix to compute the diagram.
  Diag_split_ij <- ripsDiag(
    X = dist1, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", dist="arbitrary", location = TRUE, printProgress = F)
  return(Diag_split_ij)
}
