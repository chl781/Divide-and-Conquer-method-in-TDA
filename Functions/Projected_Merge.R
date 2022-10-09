# projected Merge step



Projected_Merge <- function(Diag_split,X_split,m,gap1,gap2){
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
  
  # Allow 2 subfeatures to merge, cancelling the projected boundaries.
  
  
  bound=matrix(list(),m-1,m-1)
  bound_number=matrix(list(),m-1,m-1)
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
        bound1=c() # This indicates the boundary points
        bound2=c() # This indicates which boundaries are included.
        data1=Diag_split_ij[["cycleLocation"]][[t3[ij]]]
        IN=as.vector(data1)
        IN1=IN-lengthX # Show the indices.
        IN1=IN1[IN1>0]
        
        IN2=unique(IN1)
        gap_space=(gap3[2]-gap3[1])/lengthX*5
        lengthBound=floor(lengthX/5)+2 # Can be changed.
        
        
        # row_IN2=vector("list",length(IN2))
        # for(iRow_IN2 in 1:length(IN2)){
        #   row_IN2[[iRow_IN2]]=c(which(data1[,1]==IN2[iRow_IN2]),
        #                         which(data1[,2]==IN2[iRow_IN2]))
        # }
        
        # Need a concrete projection algorithm.
        # projection may also need a threshold.
        
        data2=X_split_ij[IN[IN<=lengthX],]
        #lengthdata=NROW(data2)
        
        # Will introduce another parameter about the gap between points.
        if(1 %in% IN1){  
          # data4=cbind(data1[data1[,1]==1+lengthX,2],data1[data1[,2]==1+lengthX,1])
          # data5=data4[data4<=lengthX]
          # data6=data4[data4>lengthX]
          # if(length(data5)>0){
          if(i!=1){
            #bound1=rbind(bound1,cbind(rep(gap3[1],sum(IN1==1)),X_split_ij[data5,2]))
            bound1=rbind(bound1,cbind(rep(gap3[1],lengthBound),seq(gap4[1],gap4[2],length.out = lengthBound)))
            bound2 = c(bound2,1)
          }else{
            #bound1=rbind(bound1,cbind(rep(gap3[2],sum(IN1==1)),X_split_ij[data5,2]))
            bound1=rbind(bound1,cbind(rep(gap3[2],lengthBound),
                                      seq(gap4[1],gap4[2],length.out = lengthBound)))
            bound2 = c(bound2,2)
          }
          #}
        }
        
        if(2 %in% IN1){ 
          # data4=cbind(data1[data1[,1]==2+lengthX,2],data1[data1[,2]==2+lengthX,1])
          # data5=data4[data4<=lengthX]
          # data6=data4[data4>lengthX]
          # if(length(data5)>0){
          if(i>1&i<m-1){
            #bound1=rbind(bound1,cbind(rep(gap3[2],sum(IN1==2)),X_split_ij[data5,2]))
            bound1=rbind(bound1,cbind(rep(gap3[2],lengthBound),
                                      seq(gap4[1],gap4[2],length.out = lengthBound)))
            bound2 = c(bound2,2)
          }else{
            if(j>1){
              #bound1=rbind(bound1,cbind(X_split_ij[data5,1],rep(gap4[1],sum(IN1==2))))
              bound1=rbind(bound1,cbind(seq(gap3[1],gap3[2],
                                            length.out = lengthBound),rep(gap4[1],
                                                                          lengthBound)))
              bound2 = c(bound2,3)
            } else{
              bound1=rbind(bound1,cbind(seq(gap3[1],gap3[2],
                                            length.out = lengthBound),rep(gap4[2],
                                                                          lengthBound)))
              bound2 = c(bound2,4)
            }
          }
          #}
        }
        
        if(3 %in% IN1){  
          # data4=cbind(data1[data1[,1]==3+lengthX,2],data1[data1[,2]==3+lengthX,1])
          # data5=data4[data4<=lengthX]
          # data6=data4[data4>lengthX]
          # if(length(data5)>0){
          if(i==1|i==m-1|j==1){
            bound1=rbind(bound1,cbind(seq(gap3[1],gap3[2],
                                          length.out = lengthBound)
                                      ,rep(gap4[2],lengthBound)))
            bound2 = c(bound2,4)
          }else{
            bound1=rbind(bound1,cbind(seq(gap3[1],gap3[2],
                                          length.out = lengthBound)
                                      ,rep(gap4[1],lengthBound)))
            bound2 = c(bound2,3)
          }
          #}
        }
        
        if(4 %in% IN1){
          # data4=cbind(data1[data1[,1]==4+lengthX,2],data1[data1[,2]==4+lengthX,1])
          # data5=data4[data4<=lengthX]
          # data6=data4[data4>lengthX]
          # if(length(data5)>0){
          bound1=rbind(bound1,cbind(seq(gap3[1],gap3[2],
                                        length.out = lengthBound)
                                    ,rep(gap4[2],lengthBound)))
          bound2 = c(bound2,4)
          #}
        }
        
        bound1 = unique(bound1)
        bound2 = unique(bound2)
        
        # adding the connection between boundaries
        # if(1%in% IN2 & 2 %in% IN2){
        #   if(i==1){
        #     if(j==1){
        #       bound1=rbind(bound1,c(gap3[2],gap4[2]))
        #     }else{
        #       bound1=rbind(bound1,c(gap3[2],gap4[1]))
        #     }
        #   }else if(i==m-1){
        #     if(j==1){
        #       bound1=rbind(bound1,c(gap3[1],gap4[2]))
        #     }else{
        #       bound1=rbind(bound1,c(gap3[1],gap4[1]))
        #     }
        #   }else{
        #     warning("Wrong Connection information")
        #   }
        # }
        # 
        # 
        # if(1%in% IN2 & 3 %in% IN2){
        #   if(i==1){
        #     if(j==1|j==m-1){
        #       warning("Wrong Connection information")
        #     } else{
        #       bound1=rbind(bound1,c(gap3[2],gap4[2]))
        #     }
        #   }else if(i==m-1){
        #     if(j==1|j==m-1){
        #       warning("Wrong Connection information")
        #     } else{
        #       bound1=rbind(bound1,c(gap3[1],gap4[2]))
        #     }
        #   }else{
        #     if(j==1){
        #       bound1=rbind(bound1,c(gap3[1],gap4[2]))
        #     }else{
        #       bound1=rbind(bound1,c(gap3[1],gap4[1]))
        #     }
        #   }
        # }
        # 
        # if(1%in% IN2 & 4 %in% IN2){
        #   bound1=rbind(bound1,c(gap3[1],gap4[2]))
        # }
        # 
        # if(2 %in% IN2 & 3%in% IN2){
        #   if(j==1){
        #     bound1=rbind(bound1,c(gap3[2],gap4[2]))
        #   }else{
        #     bound1=rbind(bound1,c(gap3[2],gap4[1]))
        #   }
        # }
        # 
        # if(2 %in% IN2 & 4 %in% IN2){
        #   bound1=rbind(bound1,c(gap3[2],gap4[2]))
        # }
        # 
        # if(3 %in% IN2 & 4 %in% IN2){
        #   warning("Wrong Connection information")
        # }
        
        # Topological recovery.
        # Naive estimate.
        
        colnames(bound1)=c("x","y")
        #data3=rbind(bound1,data2)
        DiagProj=ripsDiag(
          X = rbind(bound1,data2), maxdimension = maxdimension, maxscale = maxscale,
          library = "Dionysus", location = TRUE, printProgress = F)
        
        ind=which(DiagProj$diagram[,1]==1)
        if(sum(ind)==0){
          bound1=matrix(c(Inf,Inf,Inf,Inf),2)
          bound[[i,j]][[ij]]=bound1
        }else{
          s=which.max(DiagProj$diagram[ind,3]-DiagProj$diagram[ind,2])
          s=ind[s]
          bound1=unique(matrix(as.vector(
            DiagProj$cycleLocation[[s]]),ncol=2,byrow = F))
          Suspicious[[i,j]][[ij]]=DiagProj$cycleLocation[[s]]
          Suspicious_features1=unique(
            matrix(as.vector(DiagProj$cycleLocation[[s]]),ncol = 2))
          bound[[i,j]][[ij]]=rbind(bound1[bound1[,1]==gap3[1]|bound1[,1]==gap3[2],],
                                   bound1[bound1[,2]==gap4[1]|bound1[,2]==gap4[2],])
          bound_number[[i,j]][[ij]]=bound2
          Suspicious_features[[i,j]][[ij]]=
            Suspicious_features1[Suspicious_features1[,1]!=gap3[1]&
                                   Suspicious_features1[,1]!=gap3[2]&
                                   Suspicious_features1[,2]!=gap4[1]&
                                   Suspicious_features1[,2]!=gap4[2],]
          }
      }
      
    }
  }
  
  ## Matching
  non_empty=c()
  for(j1 in 1:(m-1)){ #j1 here
    for(i1 in 1:(m-1)){ #i1 here
      t=1
      if(length_bound_matrix[i1,j1]!=0){
        non_empty=rbind(non_empty,
                        cbind(matrix(rep(c(i1,j1),length_bound_matrix[i1,j1]),
                                     nrow = length_bound_matrix[i1,j1],
                                     byrow = T),1:length_bound_matrix[i1,j1]))
      }
    }
  }
  
  length_bound_matrix=matrix(0,m-1,m-1)
  for (i in 1 : (m-1)) {
    for (j in 1 : (m-1)) {
      length_bound_matrix[i,j]=length(bound[[i,j]])
    }
  }
  length_bound=sum(length_bound_matrix)
  
  Combined=vector("list",length(suspicious_ind))
  num=1
  for(i in 1:(length_bound-1)){
    for (j in (i+1):length_bound) {
      Ndata=floor(max(NROW(Suspicious[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]]),
                NROW(Suspicious[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]]))/5)+2
      error=max(gap1[2]-gap1[1],gap2[2]-gap2[1])/Ndata
      if(Matching1(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]]
                  ,bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]],eps = error )){
        df=Diag_split[[non_empty[i,1],non_empty[i,2]]][["cycleLocation"]]
        df1=df[[ind_suspicious[[non_empty[i,1], non_empty[i,2] ]][ non_empty[i,3] ] ]]
        df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i,1],non_empty[i,2]]])])

        Suspicious_i=X_split[[non_empty[i,1],non_empty[i,2]]][df1,]
        
        df=Diag_split[[non_empty[j,1],non_empty[j,2]]][["cycleLocation"]]
        df1=df[[ind_suspicious[[non_empty[j,1], non_empty[j,2] ]][ non_empty[j,3] ] ]]
        df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[j,1],non_empty[j,2]]])])

        Suspicious_j=X_split[[non_empty[j,1],non_empty[j,2]]][df1,]
        
        Combined[[num]]=rbind(Suspicious_i,Suspicious_j)
        print(i)
        print(j)
        num=num+1
      }
    }
  }
    
  # Allow 3 subfeatures to merge, cancelling the projected boundaries.
  if(length_bound>=3){
    for(i in 1:(length_bound-2)){
      for (j in (i+1):(length_bound-1)) {
        for(k in (j+1):length_bound)
          Ndata=floor(max(NROW(Suspicious[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]]),
                          NROW(Suspicious[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]]))/5)+2
        error=max(gap1[2]-gap1[1],gap2[2]-gap2[1])/Ndata
        list_X=list(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]],
                    bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]],
                    bound[[non_empty[k,1],non_empty[k,2]]][[non_empty[k,3] ]])
        if(Matching2(list_X,eps = error )){
          for(i1 in c(i,j,k)){
            df=Diag_split[[non_empty[i1,1],non_empty[i1,2]]][["cycleLocation"]]
            df1=df[[ind_suspicious[[non_empty[i1,1], non_empty[i1,2] ]][ non_empty[i1,3] ] ]]
            df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i1,1],non_empty[i1,2]]])])
            
            Suspicious_i=X_split[[non_empty[i1,1],non_empty[i,2]]][df1,]
            Combined[[num]]=rbind(Combined[[num]],Suspicious_i)
          }
          num=num+1
        }
      }
    }
  }
  
 return(Combined)
}