# This code is used to test how the combination behaving but cut into 4 pieces.
# Return as a matrix with full data birth time and death time and two
# cutting pieces' birth time and death times.
# n is the number of supplemental points.

DiagCir4Pieces <- function(X,gap1,gap2,maxscale,maxdimension,n){
  #Split data into 4 parts.
  X1=X[-X[,1]>gap1&-X[,2]>gap2,]
  X2=X[-X[,1]<gap1&-X[,2]>gap2,]
  X3=X[-X[,1]>gap1&-X[,2]<gap2,]
  X4=X[-X[,1]<gap1&-X[,2]<gap2,]
  
  Y1=matrix(c(rep(-gap1,n),seq(-1.1,1.1,2.2/(n-1))),nrow=n,byrow = F)
  Y2=matrix(c(seq(-1.1,1.1,2.2/(n-1)),rep(-gap2,n)),nrow=n,byrow = F)
  S=matrix(0,nrow = 9,ncol = 2)
  
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
  
  # X1
  DiagRips1 <- ripsDiag(
    X = rbind(X1,Y2[-Y2[,1]>gap1,],Y1[-Y1[,2]>gap2,]), 
    maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips1$diagram[,1]==1) ){
    index=which(DiagRips1$diagram[,1]==1)
    t1=index[which.max(DiagRips1$diagram[index,3]-DiagRips1$diagram[index,2])]
    S[2,]=c(DiagRips1$diagram[t1,2],DiagRips1$diagram[t1,3])
  }else{
    t1=-1
    S[2,]=c(0,0)
  }
  
  # X2
  DiagRips2 <- ripsDiag(
    X = rbind(X2,Y2[-Y2[,1]<gap1,],Y1[-Y1[,2]>gap2,]),
    maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips2$diagram[,1]==1) ){
    index=which(DiagRips2$diagram[,1]==1)
    t2=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
    S[3,]=c(DiagRips2$diagram[t2,2],DiagRips2$diagram[t2,3])
  }else{
    t2=-1
    S[3,]=c(0,0)
  }
  
  # X3
  DiagRips3 <- ripsDiag(
    X = rbind(X3,Y2[-Y2[,1]>gap1,],Y1[-Y1[,2]<gap2,]), 
    maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips3$diagram[,1]==1) ){
    index=which(DiagRips3$diagram[,1]==1)
    t3=index[which.max(DiagRips3$diagram[index,3]-DiagRips3$diagram[index,2])]
    S[4,]=c(DiagRips3$diagram[t3,2],DiagRips3$diagram[t3,3])
  }else{
    t3=-1
    S[4,]=c(0,0)
  }
  
  
  # X4
  DiagRips4 <- ripsDiag(
    X = rbind(X4,Y2[-Y2[,1]<gap1,],Y1[-Y1[,2]<gap2,]), 
    maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips4$diagram[,1]==1) ){
    index=which(DiagRips4$diagram[,1]==1)
    t4=index[which.max(DiagRips4$diagram[index,3]-DiagRips4$diagram[index,2])]
    S[5,]=c(DiagRips4$diagram[t4,2],DiagRips4$diagram[t4,3])
  }else{
    t4=-1
    S[5,]=c(0,0)
  }
  
  
  # Combine the data.
  if(t1!=-1){
    data1=unique(matrix(as.vector(DiagRips1[["cycleLocation"]][[t1]]),ncol=2,byrow = F))
  }else{
    data1=c()
  }
  if(t2!=-1){
    data2=unique(matrix(as.vector(DiagRips2[["cycleLocation"]][[t2]]),ncol=2,byrow = F))
  }else{
    data2=c()
  }
  if(t3!=-1){
    data3=unique(matrix(as.vector(DiagRips3[["cycleLocation"]][[t3]]),ncol=2,byrow = F))
  }else{
    data3=c()
  }
  if(t4!=-1){
    data4=unique(matrix(as.vector(DiagRips4[["cycleLocation"]][[t4]]),ncol=2,byrow = F))
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
  
  if(t1==-1|t2==-1|t3==-1|t4==-1){ # This is to see if this method could recover the feature.
    S[9,]=c(1,1)
  }
  
  return(S)
}