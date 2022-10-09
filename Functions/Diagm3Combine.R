# This function serves as projection method v3.

# Define the inputs
# 

Diagm3Combine <- function(X_split,m,Diag_split,
                          range,maxdimension,maxscale,error){
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

  # Combine the data from sub-regions by using the sub-feature data but not projected boundary.
  length_bound_matrix=matrix(0,m-1,m-1)
  for (i in 1 : (m-1)) {
    for (j in 1 : (m-1)) {
      length_bound_matrix[i,j]=length(ind_suspicious[[i,j]])
    }
  }
  length_bound=sum(length_bound_matrix)
  
  dist_bound=matrix(Inf,length_bound,length_bound)
  for (i1 in 1:length_bound ) {
      dist_bound[[i1,i1]]=0
  }
  
  
  for (i1 in 1:((m-1)^2-1) ) {
    for (j1 in (i1+1):((m-1)^2) ) {
      #for(j1 in c(i1+1,i1+m-1)) {
      if ( length_bound_matrix[i1]==0|length_bound_matrix[j1]==0 ){
        next
      }
      
      for (i2 in 1:length_bound_matrix[i1]) {
        for (j2 in 1:length_bound_matrix[j1]){
          df1=Diag_split[[i1]][["cycleLocation"]][[ind_suspicious[[i1]][i2]]]
          df1=unique(df1[as.vector(df1)<=NROW(X_split[[i1]])])
          df11=X_split[[i1]][df1,]
          if(NROW(df11)<=1|NROW(df11)!=length(df1)){
            #stop()
            warning()
          }
          df2=Diag_split[[j1]][["cycleLocation"]][[ind_suspicious[[j1]][j2]]]
          df2=unique(df2[as.vector(df2)<=NROW(X_split[[j1]])])
          df22=X_split[[j1]][df2,]
          if(NROW(df22)<=1|NROW(df22)!=length(df2)){
            #stop()
            warning()
          }
          
          if(i1>1){
            i3=sum(length_bound_matrix[1:(i1-1)])+i2
          }
          if(i1==1){
            i3=i2
          }
          j3=sum(length_bound_matrix[1:(j1-1)])+j2

          dist_bound[i3,j3]=min(dist(df11,df22))
          dist_bound[j3,i3]=dist_bound[i3,j3]
        }
      }
    }
  }
  
  # for the subfeatures within each block
  # for(i1 in 1:NROW(dist_bound) ){
  #   for (j1 in 1:NROW(dist_bound) ) {
  #     if(identical(dist_bound[i1,j1],0) & i1!=j1){
  #       dist_bound[i1,j1]=Inf
  #     }
  #   }
  # }
  
  
  for (i1 in 1:(m-1)^2 ) {
    length_bound_i1=length_bound_matrix[i1]
    if(length_bound_i1>1){
      for(i2 in 1:(length_bound_i1-1) ){
        for(j2 in (i2+1):length_bound_i1 ){
          df1=Diag_split[[i1]][["cycleLocation"]][[ind_suspicious[[i1]][i2]]]
          df1=unique(df1[as.vector(df1)<=NROW(X_split[[i1]])])
          df11=X_split[[i1]][df1,]
          if(NROW(df11)<=1|NROW(df11)!=length(df1)){
            stop()
          }
          df2=Diag_split[[i1]][["cycleLocation"]][[ind_suspicious[[i1]][j2]]]
          df2=unique(df2[as.vector(df2)<=NROW(X_split[[i1]])])
          df22=X_split[[i1]][df2,]
          if(NROW(df22)<=1|NROW(df22)!=length(df2)){
            stop()
          }
          
          if(i1>1){
            i3=sum(length_bound_matrix[1:(i1-1)])+i2
            j3=sum(length_bound_matrix[1:(i1-1)])+j2
          }
          if(i1==1){
            i3=i2
            j3=j2
          }
          dist_bound[i3,j3]=min(dist(df11,df22))
          dist_bound[j3,i3]=dist_bound[i3,j3]
        }
      }
    }
  }
  
  # Notice that the maxscale of the following function cannot be too large.
  
  diag_suspicious=ripsDiag(dist_bound,maxdimension,maxscale = maxscale,
                           library = "Dionysus",dist="arbitrary",
                           location = T)
  suspicious_ind=which(diag_suspicious$diagram[,1]==1)
  Combined=vector("list",length(suspicious_ind))
  Combined_bound=vector("list",length(suspicious_ind))
  
  # Record the index that save the data in suspicious data.
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
  
  num=1
  for(one in suspicious_ind){
    Combined_ind=unique(as.vector(diag_suspicious$cycleLocation[[one]]))
    col=1
    for(i in Combined_ind){
      df=Diag_split[[non_empty[i,1],non_empty[i,2]]][["cycleLocation"]]
      df1=df[[ind_suspicious[[non_empty[i,1], non_empty[i,2] ]][ non_empty[i,3] ] ]]
      df1=unique(df1[as.vector(df1)<=NROW(X_split[[non_empty[i,1],non_empty[i,2]]])])
      if(length(df1)<=1){
        stop()
      }
      Suspicious_i=cbind(X_split[[non_empty[i,1],non_empty[i,2]]][df1,],
                         rep(col,NROW(X_split[[non_empty[i,1],non_empty[i,2]]][df1,])))
      colnames(Suspicious_i)<-c("x","y","col")
      Combined[[num]]=rbind(Combined[[num]],Suspicious_i)
      col=col+1
    }
    num=num+1
  }
  
  # The following uses the projected method.
  Projected_Merge1 <- Projected_Merge(Diag_split,X_split,m,gap1,gap2)
  
  # Recover the estimates below.
  
  
  PD <- matrix(0,nrow=3,ncol=NROW(Combined)+NROW(Projected_Merge1))
  colnames(PD) <- c("dimension","birth","death")
  PD[,1]=1
  for (i in 1:NROW(PD) ) {
    PD[i,2]=
  }
  
  return(list(diagram=PD,Combined))
}