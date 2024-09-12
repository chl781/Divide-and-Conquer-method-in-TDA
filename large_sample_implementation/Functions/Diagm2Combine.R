# This function serves as projection method.

# Define the inputs
# 

Diagm2Combine <- function(X_split,m,Diag_split,
                          range,maxdimension,maxscale){
  ind_suspicious=matrix(list(),m-1,m-1)
  ind_subfeature=matrix(list(),m-1,m-1)
  
  
  gap1=seq(range[1,1],range[1,2],length.out =m)
  gap2=seq(range[2,1],range[2,2],length.out =m)
  
  # Extract the information of cross feature.
  
  # lapply
  # mclapply
  
  for (i in 1:(m-1)) {
    for (j in 1:(m-1)) {
      t3=c() # temporary variable for ind_suspicious
      t4=c() # temporary variable for ind_subfeature
      lengthX=NROW(X_split[[i,j]])
      Diag_split_ij=Diag_split[[i,j]]
      index=which(Diag_split_ij$diagram[,1]==1)
      if(length(index)==0){
        next
      }
      for (ij in 1:length(index)){
        if((lengthX+1) %in% as.vector(Diag_split_ij[["cycleLocation"]][[index[ij]]])|
           (lengthX+2) %in% as.vector(Diag_split_ij[["cycleLocation"]][[index[ij]]])|
           (lengthX+3) %in% as.vector(Diag_split_ij[["cycleLocation"]][[index[ij]]])|
           (lengthX+4) %in% as.vector(Diag_split_ij[["cycleLocation"]][[index[ij]]])){
          ###################
          # skip ij if birth time is 0.
          if(Diag_split_ij$diagram[index[ij],2]!=0){
            t3=c(t3,index[ij])
          }
          ###################
        }else{
          t4=c(t4,index[ij])
        }
      }
      if(!is.null(t3)){
        ind_suspicious[[i,j]]=t3
      }
      if(!is.null(t4)){
        ind_subfeature[[i,j]]=t4
      }
    }
  }
  
  # Project information.
  bound=matrix(list(),m-1,m-1)
  Suspicious=matrix(list(),m-1,m-1)
  Suspicious_features=matrix(list(),m-1,m-1)
  for (i in 1:(m-1)) {
    for (j in 1:(m-1)) {
      if(length(ind_suspicious[[i,j]])==0){
        next
      }
      X_split_ij=X_split[[i,j]]
      #colnames(X_split_ij)=NULL
      Diag_split_ij=Diag_split[[i,j]]
      gap3=c(gap1[i],gap1[i+1])
      gap4=c(gap2[j],gap2[j+1])
      
      t3=ind_suspicious[[i,j]]
      
      bound[[i,j]]=matrix(list(),1,length(t3))
      Suspicious[[i,j]]=matrix(list(),1,length(t3))
      Suspicious_features[[i,j]]=matrix(list(),1,length(t3))
      # This discussion is for general.
      lengthX=NROW(X_split[[i,j]])
      
      for(ij in 1:length(t3)){
        bound1=c()
        IN=unique(as.vector(Diag_split_ij[["cycleLocation"]][[t3[ij]]]))
        IN1=IN-lengthX # Show the indices.
        IN1=IN1[IN1>0]
        
        # Need a concrete projection algorithm.
        # projection may also need a threshold.
        data2=X_split_ij[IN[IN<=lengthX],]
        lengthdata=NROW(data2)
        if(1 %in% IN1){  
          if(i!=1){
            bound1=rbind(bound1,cbind(rep(gap3[1],lengthdata),data2[,2]))
          }else{
            bound1=rbind(bound1,cbind(rep(gap3[2],lengthdata),data2[,2]))
          }
        }
        if(2 %in% IN1){  
          if(i>1&i<m-1){
            bound1=rbind(bound1,cbind(rep(gap3[2],lengthdata),data2[,2]))
          }else{
            if(j>1){
              bound1=rbind(bound1,cbind(data2[,1],rep(gap4[1],lengthdata)))
            } else{
              bound1=rbind(bound1,cbind(data2[,1],rep(gap4[2],lengthdata)))
            }
          }
        }
        
        if(3 %in% IN1){  
          if(i==1|i==m-1||j==1){
            bound1=rbind(bound1,cbind(data2[,1],rep(gap4[2],lengthdata)))
          }else{
            bound1=rbind(bound1,cbind(data2[,1],rep(gap4[1],lengthdata)))
          }
        }
        
        if(4 %in% IN1){  
          bound1=rbind(bound1,cbind(data2[,1],rep(gap4[2],lengthdata)))
        }
        
        # Topological recovery.
        # Naive estimate.
        colnames(bound1)=c("x","y")
        data3=rbind(bound1,data2)
        DiagProj=ripsDiag(
          X = data3, maxdimension = maxdimension, maxscale = maxscale,
          library = "Dionysus", location = TRUE, printProgress = F)
        
        ind=which(DiagProj$diagram[,1]==1)
        s=which.max(DiagProj$diagram[ind,3]-DiagProj$diagram[ind,2])
        s=ind[s]
        bound1=unique(matrix(as.vector(
          DiagProj$cycleLocation[[s]]),ncol=2,byrow = F))
        bound[[i,j]][[ij]]=rbind(bound1[bound1[,1]==gap3[1]|bound1[,1]==gap3[2],],
                                 bound1[bound1[,2]==gap4[1]|bound1[,2]==gap4[2],])
        Suspicious[[i,j]][[ij]]=DiagProj$cycleLocation[[s]]
        Suspicious_features1=unique(
          matrix(as.vector(DiagProj$cycleLocation[[s]]),ncol = 2))
        Suspicious_features[[i,j]][[ij]]=
          Suspicious_features1[Suspicious_features1[,1]!=gap3[1]&
                                 Suspicious_features1[,1]!=gap3[2]&
                                 Suspicious_features1[,2]!=gap4[1]&
                                 Suspicious_features1[,2]!=gap4[2],]
      }
      
    }
  }
  

  
  ################# Matching
  
  #Suspicious_Combined=list()
  #Suspicious_Index=list()
  
  #Suspicious_num=1
  
  # for (j1 in 1:(m-1)) {
  #   for (i1 in 1:(m-1)) {
  #     # Propogate from left to right.
  #     bound_ij=bound[[i1,j1]]
  #     if(length(bound_ij)==0){
  #       next
  #     }
  #     
  #     for(k in 1:length(bound_ij)){
  #       bound_ij_k=bound_ij[[k]]
  #       if(i1 != m){
  #         #Match bound_ij[[k]]
  #         if(gap1[i1+1] %in% bound_ij_k[,1]){
  #           bound_k=bound_ij_k[bound_ij_k[,1]==gap1[i1+1],]
  #           bound_i1j=bound[[i1+1,j1]]
  #           if(length(bound_i1j)>=1){
  #             for(l in 1:length(bound_i1j)){
  #               bound_i1j_l=bound_i1j[[l]]
  #               bound_l=bound_i1j_l[bound_i1j_l[,1]==gap1[i1+1],]
  #               # Check if bound_l and bound_k can match.
  #               label=Matching(bound_l,bound_k,eps)
  #               if(label){
  #                 Suspicious_Combined[[Suspicious_num]]=
  #                 Suspicious_Index[[Suspicious_num]]=Suspicious_Index[[Suspicious_num]] 
  #                 Suspicious_num=Suspicious_num+1
  #                 
  #               }
  #             }
  #           }
  #         }
  #       }
  #       
  #       if(j1 != m){
  #         
  #       }
  #       
  #     }
  #   }
  # }
  ################
  
  # Alternative combination method.
  # Designed for large and middle features.
  # Construct the distance matrix.
  
  length_bound_matrix=matrix(0,m-1,m-1)
  for (i in 1 : (m-1)) {
    for (j in 1 : (m-1)) {
      length_bound_matrix[i,j]=length(bound[[i,j]])
    }
  }
  length_bound=sum(length_bound_matrix)
  dist_bound=matrix(0,length_bound,length_bound)
  
  for (i1 in 1:((m-1)^2-1) ) {
    for (j1 in (i1+1):(m-1)^2) {
    #for(j1 in c(i1+1,i1+m-1)) {
      if (length_bound_matrix[i1]==0|length_bound_matrix[j1]==0){
        next
      }
      #################### Can only accept combining two nearby regions.
      #if (i1 %% (m-1)==0){
      #  if(j1==i1+1){
      #    next
      #  }
      #}
      #if(j1>(m-1)^2){
      #  next
      #}
      ####################
      for (i2 in 1:length_bound_matrix[i1]) {
        for (j2 in 1:length_bound_matrix[j1]){
          df1=bound[[i1]][[i2]]
          df2=bound[[j1]][[j2]]
          if(i1>1){
            i3=sum(length_bound_matrix[1:(i1-1)])+i2
          }
          if(i1==1){
            i3=i2
          }
          j3=sum(length_bound_matrix[1:(j1-1)])+j2
          dist_bound[i3,j3]=min(dist(df1,df2))
          dist_bound[j3,i3]=dist_bound[i3,j3]
        }
      }
    }
  }
  
  diag_suspicious=ripsDiag(dist_bound,maxdimension,maxscale,
                           library = "Dionysus",dist="arbitrary",
                           location = T)
  suspicious_ind=which(diag_suspicious$diagram[,1]==1)
  Combined=vector("list",length(suspicious_ind))
  
  # Record the index that save the data in suspicious data.
  non_empty=c()
  for(i1 in 1:(m-1)){
    for(j1 in 1:(m-1)){
      if(length(Suspicious_features[[i1,j1]])!=0){
        non_empty=rbind(non_empty,
                        matrix(rep(c(i1,j1),length(Suspicious_features[[i1,j1]])),
                               nrow = length(Suspicious_features[[i1,j1]]),byrow = T))
      }
    }
  }

  num=1
  for(one in suspicious_ind){
    Combined_ind=unique(as.vector(diag_suspicious$cycleLocation[[one]]))
    for(i in Combined_ind){
      Suspicious_i=Suspicious_features[[non_empty[i,1],non_empty[i,2]]][[1]]
      Combined[[num]]=rbind(Combined[[num]],Suspicious_i)
    }
    num=num+1
  }
  return(Combined)
}