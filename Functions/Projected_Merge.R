# projected Merge step
# Check if some suspicious features can cancel the projected boundaries
# Will ignore if the length of boundary is smaller than 3. (too insignificant)


Projected_Merge <- function(ind_suspicious,bound,range,Diag_split,X_split,Suspicious,non_empty,
                            m,gap1,gap2){
  # Allow 2 subfeatures to merge, cancelling the projected boundaries.
  
  ## Matching
  
  length_bound_matrix=matrix(0,m-1,m-1)
  for (i in 1 : (m-1)) {
    for (j in 1 : (m-1)) {
      length_bound_matrix[i,j]=length(bound[[i,j]])
    }
  }
  length_bound=sum(length_bound_matrix)
  
  Combined=vector("list")
  Combined_diag_indices=vector("list")
  num=1
  for(i in 1:(length_bound-1)){
    if(NROW(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]])<=4){
      next
    }
    for (j in (i+1):length_bound) {
      if(NROW(bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]])<=4){
        next
      }
      Ndata=max(floor(max(NROW(Suspicious[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]]),
                NROW(Suspicious[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]]))/5)+2)
      error=max(gap1[2]-gap1[1],gap2[2]-gap2[1])/Ndata
      #if(Matching1(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]]
       #           ,bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]],eps = error )){
      list_X=list(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]],
                  bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]])
      
      if(Matching2(list_X,eps = error )){
        df=Diag_split[[non_empty[i,1],non_empty[i,2]]][["cycleLocation"]]
        df1=df[[ind_suspicious[[non_empty[i,1], non_empty[i,2] ]][ non_empty[i,3] ] ]]
        df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i,1],non_empty[i,2]]])])

        Suspicious_i=X_split[[non_empty[i,1],non_empty[i,2]]][df1,]
        
        df=Diag_split[[non_empty[j,1],non_empty[j,2]]][["cycleLocation"]]
        df1=df[[ind_suspicious[[non_empty[j,1], non_empty[j,2] ]][ non_empty[j,3] ] ]]
        df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[j,1],non_empty[j,2]]])])

        Suspicious_j=X_split[[non_empty[j,1],non_empty[j,2]]][df1,]
        
        Combined[[num]]=rbind(Suspicious_i,Suspicious_j)
        Combined_diag_indices=append(Combined_diag_indices,list(rbind(non_empty[j,],non_empty[i,]) ) )
        num=num+1
      }
    }
  }
    
  # Allow 3 subfeatures to merge, cancelling the projected boundaries.
  if(length_bound>=3){
    # To prevent the trivial case.
    for(i in 1:(length_bound-2)){
      if(NROW(bound[[non_empty[i,1],non_empty[i,2]]][[non_empty[i,3] ]])<=4){
        next
      }
      for (j in (i+1):(length_bound-1)) {
        if(NROW(bound[[non_empty[j,1],non_empty[j,2]]][[non_empty[j,3] ]])<=4){
          next
        }
        for(k in (j+1):length_bound){
          if(NROW(bound[[non_empty[k,1],non_empty[k,2]]][[non_empty[k,3] ]])<=4){
            next
          }
        
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
              Combined_diag_indices[[num]]=rbind(Combined_diag_indices[[num]],
                                               non_empty[j,],non_empty[i,],non_empty[k,])
            }
            num=num+1
          }
        }
      }
    }
  }
  
 return(list(Combined,Combined_diag_indices))
}