# Inputs:
# X_split is the data separate into several subregions.
# m is the number of subregions to divide.
# Diag_split is the persistent diagram with cycles into different subregions.
# range is the range of data that we consider.
# maxdimension and maxscale are inputs to ripsDiag
# error is the error bound allowed to use in the cancellation method.
# if_retrieve_suspicious_data is a boolean flag indicating whether to retrieve the suspicious data.

Diagm3Combine <- function(X_split,m,Diag_split,
                          range,maxdimension,maxscale,error,if_retrieve_suspicious_data=TRUE){
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
  # Build X_suspicious and its source-index map at the same time.  Previously,
  # X_suspicious was filtered by drop_list(), whereas non_empty was built from
  # the unfiltered length_bound_matrix.  As soon as one short cycle was removed,
  # row k of non_empty no longer described X_suspicious[[k]].
  X_suspicious <- list()
  non_empty_rows <- list()
  non_empty_all_rows <- list()

  # Outer j / inner i preserves the column-major order used when the matrices
  # above were created from mapply()/sapply().
  for (j1 in seq_len(m - 1L)) {
    for (i1 in seq_len(m - 1L)) {
      local_diag_ids <- ind_suspicious[[i1, j1]]

      if (length(local_diag_ids) != length_bound_matrix[i1, j1]) {
        stop(sprintf(
          "Inconsistent suspicious-feature count in block (%d, %d): ind_suspicious has %d but length_bound_matrix records %d.",
          i1, j1, length(local_diag_ids), length_bound_matrix[i1, j1]
        ))
      }

      for (local_feature in seq_along(local_diag_ids)) {
        all_id <- length(non_empty_all_rows) + 1L
        non_empty_all_rows[[all_id]] <- c(
          block_i = i1,
          block_j = j1,
          local_feature = local_feature
        )

        diag_id <- local_diag_ids[[local_feature]]
        point_ids <- unique(as.integer(as.vector(
          Diag_split[[i1, j1]][["cycleLocation"]][[diag_id]]
        )))
        point_ids <- point_ids[
          !is.na(point_ids) &
            point_ids >= 1L &
            point_ids <= NROW(X_split[[i1, j1]])
        ]

        candidate_cycle <-
          X_split[[i1, j1]][point_ids, , drop = FALSE]
        candidate_cycle <- drop_list(candidate_cycle)

        # Apply the original filter, but discard the source-index row together
        # with the cycle so that the two objects cannot become misaligned.
        if (is.null(candidate_cycle)) {
          next
        }

        new_id <- length(X_suspicious) + 1L
        X_suspicious[[new_id]] <- candidate_cycle
        non_empty_rows[[new_id]] <- c(
          block_i = i1,
          block_j = j1,
          local_feature = local_feature
        )
      }
    }
  }

  if (length(non_empty_rows) > 0L) {
    non_empty <- do.call(rbind, non_empty_rows)
    storage.mode(non_empty) <- "integer"
  } else {
    non_empty <- matrix(
      integer(0), nrow = 0L, ncol = 3L,
      dimnames = list(NULL, c("block_i", "block_j", "local_feature"))
    )
  }

  # Projected_Merge() historically receives the complete map, including short
  # cycles. Keep that contract separate from the filtered Rips-vertex map.
  if (length(non_empty_all_rows) > 0L) {
    non_empty_all <- do.call(rbind, non_empty_all_rows)
    storage.mode(non_empty_all) <- "integer"
  } else {
    non_empty_all <- matrix(
      integer(0), nrow = 0L, ncol = 3L,
      dimnames = list(NULL, c("block_i", "block_j", "local_feature"))
    )
  }

  stopifnot(length(X_suspicious) == NROW(non_empty))
  stopifnot(NROW(non_empty_all) == length_bound)
  
  ##
  ## dist construction
  length_dist_bound = length(X_suspicious)
  if (length_dist_bound > 0L) {
    dist_bound = matrix(
      mapply(
        FUN = dist_construct,
        X = rep(X_suspicious, times = length_dist_bound),
        Y = rep(X_suspicious, each = length_dist_bound)
      ),
      nrow = length_dist_bound,
      ncol = length_dist_bound
    )
  } else {
    dist_bound = matrix(numeric(0), nrow = 0L, ncol = 0L)
  }
  
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

    
  # Notice that the maxscale in the following function cannot be too large.
  # Running the Rips filtration based on the dist matrix among different suspicious features.
  
    if (length_dist_bound > 0L) {
      diag_suspicious = ripsDiag(
        dist_bound, maxdimension, maxscale = maxscale,
        library = "Dionysus", dist = "arbitrary",
        location = TRUE
      )
      suspicious_ind = which(diag_suspicious$diagram[, 1] == 1)
    } else {
      # ripsDiag() cannot operate on a 0 x 0 distance matrix.
      diag_suspicious = list(
        diagram = matrix(numeric(0), nrow = 0L, ncol = 3L),
        cycleLocation = list()
      )
      suspicious_ind = integer(0)
    }
    Combined = vector("list",length(suspicious_ind))
    Combined_diag_indices = vector("list",length(suspicious_ind))
    Combined_bound = vector("list",length(suspicious_ind))

    ##### Retrieve the data

    #The below is to retrieve the data splited for distance method.
    if (if_retrieve_suspicious_data){
    num=1
    for(one in suspicious_ind){
      Combined_ind = unique(as.vector(diag_suspicious$cycleLocation[[one]]))
      if (anyNA(Combined_ind) ||
          any(Combined_ind < 1L | Combined_ind > NROW(non_empty))) {
        stop(sprintf(
          "ripsDiag returned a cycle vertex outside non_empty: valid range is 1:%d; observed range is %s:%s.",
          NROW(non_empty), min(Combined_ind, na.rm = TRUE),
          max(Combined_ind, na.rm = TRUE)
        ))
      }
      Combined_ind <- as.integer(Combined_ind)
      #col = 1
      for(i in Combined_ind){
        df = Diag_split[[non_empty[i,1],non_empty[i,2]]][["cycleLocation"]] # Retrieve the index
        df1 = df[[ind_suspicious[[non_empty[i,1], non_empty[i,2] ]][ non_empty[i,3] ] ]] 
        df1 = unique(as.integer(as.vector(df1)))
        df1 = df1[
          !is.na(df1) & df1 >= 1L &
            df1 <= NROW(X_split[[non_empty[i,1],non_empty[i,2]]])
        ]
        if(length(df1) <= 1){
          stop(sprintf(
            "Mapped suspicious cycle %d in block (%d, %d) has fewer than two valid points.",
            non_empty[i,3], non_empty[i,1], non_empty[i,2]
          ))
        }
        Suspicious_i = cbind(X_split[[non_empty[i,1],non_empty[i,2]]][df1, , drop = FALSE]
                            #,rep(col,NROW(X_split[[non_empty[i,1],non_empty[i,2]]][df1,]))
                            ) # Retrieve the data
        
        #colnames(Suspicious_i) <- c("x","y","col")
        colnames(Suspicious_i) <- c("x","y")
        Combined[[num]] = rbind(Combined[[num]],Suspicious_i)
        Combined_diag_indices[[num]] = rbind(
          Combined_diag_indices[[num]],
          non_empty[i, , drop = FALSE]
        ) # This saves all of points constructing the loops.
        #col=col+1
      }
      num = num+1
    }
  }
  ##############  
    
  ##### The following code tries to use cancellation method to merge. Only allow using
  ##### 2 or 3 suspicious features to merge.
  
  # Create bound for suspicious features.
  Suspicious_bound <- boundFind(
    ind_suspicious = ind_suspicious,
    Diag_split = Diag_split,
    X_split = X_split,
    m = m,
    gap1 = gap1,
    gap2 = gap2,
    eps = error,
    maxdimension = maxdimension,
    maxscale = maxscale
  )
  bound <- Suspicious_bound$bound
  Suspicious <- Suspicious_bound$Suspicious
  
  # The following uses the projected method.
  Projected_Merge_ <- Projected_Merge(ind_suspicious,bound,range,Diag_split,X_split,Suspicious,non_empty_all
                                      ,m,gap1,gap2)
  Projected_Merge1=Projected_Merge_[[1]]
  
  
  # Combine the lists
  Combined1 <- append(Combined,Projected_Merge1)
  Combined_diag_indices2=Projected_Merge_[[2]]
  Combined_diag_indices1=append(Combined_diag_indices,Combined_diag_indices2)
  t<-integer(0)
  for(i in seq_along(Combined1)){
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
  for (i in seq_len(NROW(PD)) ) {
    PD[i,2] = BirthRecal2(Diag_split,Combined_diag_indices1[[i]],ind_suspicious)
    PD[i,3] = DeathRecal0(unique(Combined1[[i]])[,c(1,2)]) # This is a specific death estimate.
  }
  
  # Combine the combined PD and PD for sub-features.
  PD_=rbind(subfeature_PD,PD)
  # The following is the representative points.
  cycle_points=append(subfeature_cycles,Combined1)
  
  return(list(diagram=PD_,cycle_points))
}
