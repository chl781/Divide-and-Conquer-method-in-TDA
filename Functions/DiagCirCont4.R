# This code is to provide a cut into 4 examples, but with supplemental lines but not points.

DiagCont4 <- function(X,gap1,gap2,maxscale,maxdimension){
  #Split data into 4 parts.
  X1=X[-X[,1]>gap1&-X[,2]>gap2,]
  X2=X[-X[,1]<gap1&-X[,2]>gap2,]
  X3=X[-X[,1]>gap1&-X[,2]<gap2,]
  X4=X[-X[,1]<gap1&-X[,2]<gap2,]
  
  S=matrix(0,nrow = 8,ncol = 2)
  
  # Full data
  DiagRips <- ripsDiag(
    X = X, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips$diagram[,1]==1) ){
    index=which(DiagRips$diagram[,1]==1)
    t=index[which.max(DiagRips$diagram[index,3]-DiagRips$diagram[index,2])]
    S[1,]=c(DiagRips$diagram[t,2],DiagRips$diagram[t,3])
  }else{
    S[1,]=c(0,0)
  }
  
  ####### Begin to split the data into parts.
  
  # X1
  ## Data generation
  n1=nrow(X1)
  dist1=matrix(0,n1+2,n1+2)
  for (i in 1:n1) {
    for(j in i:n1){
      dist1[i,j]=dist1[j,i]=norm(X1[i,]-X1[j,],type="2")
    }
  }
  i=n1+1
  for(j in 1:n1){
    dist1[i,j]=dist1[j,i]=abs(X1[j,1])
  }
  
  i=n1+2
  for(j in 1:n1){
    dist1[i,j]=dist1[j,i]=abs(X1[j,2])
  }
  
  
  # X1 # Change X to be distance matrix
  DiagRips1 <- ripsDiag(
    X = dist1, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus",dist="arbitrary", location = TRUE, printProgress = F)
  if( any(DiagRips1$diagram[,1]==1) ){
    index=which(DiagRips1$diagram[,1]==1)
    t1=index[which.max(DiagRips1$diagram[index,3]-DiagRips1$diagram[index,2])]
    S[2,]=c(DiagRips1$diagram[t1,2],DiagRips1$diagram[t1,3])
  }else{
    t1=-1
    S[2,]=c(0,0)
  }
  
  #######
  
  # X2
  n2=nrow(X2)
  dist2=matrix(0,n2+2,n2+2)
  for (i in 1:n2) {
    for(j in i:n2){
      dist2[i,j]=dist2[j,i]=norm(X2[i,]-X2[j,],type="2")
    }
  }
  i=n2+1
  for(j in 1:n2){
    dist2[i,j]=dist2[j,i]=abs(X2[j,1])
  }
  
  i=n2+2
  for(j in 1:n2){
    dist2[i,j]=dist2[j,i]=abs(X2[j,2])
  }
  
  # X2
  DiagRips2 <- ripsDiag(
    X = dist2, maxdimension = maxdimension, dist = "arbitrary", maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips2$diagram[,1]==1) ){
    index=which(DiagRips2$diagram[,1]==1)
    t2=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
    S[3,]=c(DiagRips2$diagram[t2,2],DiagRips2$diagram[t2,3])
  }else{
    t2=-1
    S[3,]=c(0,0)
  }
  
  
  #########
  
  # X3
  n3=nrow(X3)
  dist3=matrix(0,n3+2,n3+2)
  for (i in 1:n3) {
    for(j in i:n3){
      dist3[i,j]=dist3[j,i]=norm(X3[i,]-X3[j,],type="2")
    }
  }
  i=n3+1
  for(j in 1:n3){
    dist3[i,j]=dist3[j,i]=abs(X3[j,1])
  }
  
  i=n3+2
  for(j in 1:n3){
    dist3[i,j]=dist3[j,i]=abs(X3[j,2])
  }
  
  # X3
  DiagRips3 <- ripsDiag(
    X = dist3, maxdimension = maxdimension, dist = "arbitrary", maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips3$diagram[,1]==1) ){
    index=which(DiagRips3$diagram[,1]==1)
    t3=index[which.max(DiagRips3$diagram[index,2]-DiagRips3$diagram[index,2])]
    S[4,]=c(DiagRips3$diagram[t3,2],DiagRips3$diagram[t3,3])
  }else{
    t3=-1
    S[4,]=c(0,0)
  }
  
  #######
  
  
  # X4
  n4=nrow(X4)
  dist4=matrix(0,n4+2,n4+2)
  for (i in 1:n4) {
    for(j in i:n4){
      dist4[i,j]=dist4[j,i]=norm(X4[i,]-X4[j,],type="2")
    }
  }
  i=n4+1
  for(j in 1:n4){
    dist4[i,j]=dist4[j,i]=abs(X4[j,1])
  }
  
  i=n4+2
  for(j in 1:n4){
    dist4[i,j]=dist4[j,i]=abs(X4[j,2])
  }
  
  # X4
  DiagRips4 <- ripsDiag(
    X = dist4, maxdimension = maxdimension, dist = "arbitrary", maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips4$diagram[,1]==1) ){
    index=which(DiagRips4$diagram[,1]==1)
    t4=index[which.max(DiagRips4$diagram[index,2]-DiagRips4$diagram[index,2])]
    S[5,]=c(DiagRips4$diagram[t4,2],DiagRips4$diagram[t4,3])
  }else{
    t4=-1
    S[5,]=c(0,0)
  }
  
  #####################
  # Combine the data.
  if(t1!=-1){
    X1=rbind(X1,c(-gap1,-gap2),c(-gap1,-gap2))
    data1=unique(X1[as.vector(DiagRips1[["cycleLocation"]][[t1]]),])
  }else{
    data1=c()
  }
  if(t2!=-1){
    X2=rbind(X2,c(-gap1,-gap2),c(-gap1,-gap2)) # Supplemental points. Not used
    data2=unique(X2[as.vector(DiagRips2[["cycleLocation"]][[t2]]),])
  }else{
    data2=c()
  }
  if(t3!=-1){
    X3=rbind(X3,c(-gap1,-gap2),c(-gap1,-gap2)) # Supplemental points. Not used
    data3=unique(X3[as.vector(DiagRips3[["cycleLocation"]][[t3]]),])
  }else{
    data3=c()
  }
  if(t4!=-1){
    X4=rbind(X4,c(-gap1,-gap2),c(-gap1,-gap2)) # Supplemental points. Not used
    data4=unique(X4[as.vector(DiagRips4[["cycleLocation"]][[t4]]),])
  }else{
    data4=c()
  }
  
  # Here is a problem, because if there is no complete half of the data, then we cannot
  # use the sophisticated method to produce it.
  
  # Notice that this procedure do not check for the coincidence for the representative points.
  data=rbind(data1,data2,data3,data4)
  data=data[data[,1]!=-gap1&data[,2]!=-gap2,]
  
  # Based on the new data, we analyze the estimated birth time and death time.
  # Re-estimate.
  DiagRips5 <- ripsDiag(
    X = data, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips5$diagram[,1]==1) ){
    index=which(DiagRips5$diagram[,1]==1)
    t5=index[which.max(DiagRips5$diagram[index,3]-DiagRips5$diagram[index,2])]
    S[6,]=c(DiagRips5$diagram[t5,2],DiagRips5$diagram[t5,3])
  }else{
    S[6,]=c(0,0)
  }
  
  # Added estimates.
  # Notice that I have changed the birth estimates to be simpler.
  
  
  BirthRecal=max(S[2,1],S[3,1],S[4,1],S[5,1]) 
  # This birth estimate is for general.
  # Should do it in a more sophisticated way.
  Deathrecal=DeathRecal_Circle(data,1,2,3)
  S[7,]=c(BirthRecal,Deathrecal)
  
  # Approximate estimate based on whole data
  gmra = gmra.create.ipca(X, eps=0, dim=1, maxKids=1, stop=4)
  res <- multiscale.rips(gmra, maxD = 1)
  
  t2<-which.max(res$diagram[res$diagram[,4]==1,2]- 
                  res$diagram[res$diagram[,4]==1,1])
  if(sum(res$diagram[res$diagram[,4]==1,])==0){
    MSE_Linear_birth=0
    MSE_Linear_death=0
  }else{
    MSE_Linear_birth=res$diagram[which(res$diagram[,4]==1)[t2],1]
    MSE_Linear_death=res$diagram[which(res$diagram[,4]==1)[t2],2]
  }
  
  S[8,]=c(MSE_Linear_birth,MSE_Linear_death)
  
  return(S)
}