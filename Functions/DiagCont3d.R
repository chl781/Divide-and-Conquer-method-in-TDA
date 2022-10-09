# This code is for shperical structure be recovered.
# Especially, adding continuous structure.


DiagCont3d <- function(X,gap,maxscale,maxdimension){
  X1=X[-X[,1]>gap,]
  X2=X[-X[,1]<gap,]
 
  S=matrix(0,nrow = 6,ncol = 2)
  # Full data
  DiagRips <- ripsDiag(
    X = X, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips$diagram[,1]==2) ){ # dimension goes to 2.
    index=which(DiagRips$diagram[,1]==2)
    t=index[which.max(DiagRips$diagram[index,3]-DiagRips$diagram[index,2])]
    S[1,]=c(DiagRips$diagram[t,2],DiagRips$diagram[t,3])
  }else{
    S[1,]=c(0,0)
  }
  
  # X1
  n1=nrow(X1)
  dist1=matrix(0,n1+1,n1+1)
  for (i in 1:n1) {
    for(j in i:n1){
      dist1[i,j]=dist1[j,i]=norm(X1[i,]-X1[j,],type="2")
    }
  }
  
  i=n1+1
  for(j in 1:n1){
    dist1[i,j]=dist1[j,i]=abs(X1[j,1])
  }
  
  #i=n1+2
  #for(j in 1:n1){
  #  dist1[i,j]=dist1[j,i]=abs(X1[j,2])
  #}
  
  # X1 # Change X to be distance matrix
  DiagRips1 <- ripsDiag(
    X = dist1, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus",dist="arbitrary", location = TRUE, printProgress = F)
  if( any(DiagRips1$diagram[,1]==1) ){
    index=which(DiagRips1$diagram[,1]==2)
    t1=index[which.max(DiagRips1$diagram[index,3]-DiagRips1$diagram[index,2])]
    S[2,]=c(DiagRips1$diagram[t1,2],DiagRips1$diagram[t1,3])
  }else{
    t1=-1
    S[2,]=c(0,0)
  }
  
  # X2
  
  n2=nrow(X2)
  dist2=matrix(0,n2+1,n2+1)
  for (i in 1:n2) {
    for(j in i:n2){
      dist2[i,j]=dist2[j,i]=norm(X2[i,]-X2[j,],type="2")
    }
  }
  i=n2+1
  for(j in 1:n2){
    dist2[i,j]=dist2[j,i]=abs(X2[j,1])
  }
  
  #i=n2+2
  #for(j in 1:n2){
  #  dist2[i,j]=dist2[j,i]=abs(X2[j,2])
  #}
  
  # X2
  DiagRips2 <- ripsDiag(
    X = dist2, maxdimension = maxdimension, dist = "arbitrary", maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips2$diagram[,1]==1) ){
    index=which(DiagRips2$diagram[,1]==2)# Change it to dimension = 2.
    t2=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
    S[3,]=c(DiagRips2$diagram[t2,2],DiagRips2$diagram[t2,3])
  }else{
    t2=-1
    S[3,]=c(0,0)
  }
  
  
  # Combine the data.
  if(t1!=-1){
    X1=rbind(X1,c(-gap,-gap,-gap))
    data1=unique(X1[as.vector(DiagRips1[["cycleLocation"]][[t1]]),])
  }else{
    data1=c()
  }
  if(t2!=-1){
    X2=rbind(X2,c(-gap,-gap,-gap)) # Supplemental points. Not used
    data2=unique(X2[as.vector(DiagRips2[["cycleLocation"]][[t2]]),])
  }else{
    data2=c()
  }
  
  # Here is a problem, because if there is no complete half of the data, then we cannot
  # use the sophisticated method to produce it.
  
  # There is another potential issue that the data does not necessarily have the 
  # same representative points. Otherwise, it is just a test program.
  
  # if(all(data1[data1[1,]==-gap,],data2[data2[1,]==-gap,]))
  
  data=rbind(data1,data2)
  data=data[data[,1]!=-gap,]
  
  # Based on the new data, we analyze the estimated birth time and death time.
  DiagRips3 <- ripsDiag(
    X = data, maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips3$diagram[,1]==2) ){
    index=which(DiagRips3$diagram[,1]==2)
    t3=index[which.max(DiagRips3$diagram[index,3]-DiagRips3$diagram[index,2])]
    S[4,]=c(DiagRips3$diagram[t3,2],DiagRips3$diagram[t3,3])
  }else{
    S[4,]=c(0,0)
  }
  
  ########## The above are about re-estimate and the following is about using approximation and 
  ########## Easy estimate method to estimate it.
  
  ## Easy estimate
  
  #BirthRecal=max(BirthRecal(data1,data2),S[2,1],S[3,1]) # This birth estimate is for general.
  BirthRecal=max(S[2,1],S[3,1])
  seed=c(1,2,3,4) # Can be changed to random numbers.
  # seed=ceiling(seed)
  Deathrecal=DeathRecal_Sphere(data[seed[1],],data[seed[2],],
                               data[seed[3],],data[seed[4],])
  S[5,]=c(BirthRecal,Deathrecal)
  
  # Approximate estimate based on whole data
  gmra = gmra.create.ipca(X, eps=0, dim=maxdimension, maxKids=1, stop=4)
  res <- multiscale.rips(gmra, maxD = maxdimension)
  
  t2<-which.max(res$diagram[res$diagram[,4]==2,2]- #Change dimension to 2
                  res$diagram[res$diagram[,4]==2,1])
  if(sum(res$diagram[res$diagram[,4]==2,])==0){
    MSE_Linear_birth=0
    MSE_Linear_death=0
  }else{
    MSE_Linear_birth=res$diagram[which(res$diagram[,4]==2)[t2],1]
    MSE_Linear_death=res$diagram[which(res$diagram[,4]==2)[t2],2]
  }
  
  S[6,]=c(MSE_Linear_birth,MSE_Linear_death)
  
  return(S)
}
