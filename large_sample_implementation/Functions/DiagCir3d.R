# 3D Data and supplement data as a plane
# This code is used to test how the combination behaving.
# Return as a matrix with full data birth time and death time and two
# cutting pieces' birth time and death times.
# n is the number of supplemental points.

# Notice that we can use the sphere information to do the easy calculation.

# This program is only used for cutting in half.

DiagCir3d <- function(X,gap,maxscale,maxdimension,n){
  X1=X[-X[,1]>gap,]
  X2=X[-X[,1]<gap,]
  Y=matrix(c(rep(-gap,n^2),rep(seq(-1.1,1.1,2.2/(n-1)),n),
             rep(seq(-1.1,1.1,2.2/(n-1)),each=n)),
           ncol=3,byrow = F)
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
  DiagRips1 <- ripsDiag(
    X = rbind(X1,Y), maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips1$diagram[,1]==2) ){
    index=which(DiagRips1$diagram[,1]==2)
    t1=index[which.max(DiagRips1$diagram[index,3]-DiagRips1$diagram[index,2])]
    S[2,]=c(DiagRips1$diagram[t1,2],DiagRips1$diagram[t1,3])
  }else{
    t1=-1
    S[2,]=c(0,0)
  }
  
  # X2
  DiagRips2 <- ripsDiag(
    X = rbind(X2,Y), maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips2$diagram[,1]==2) ){
    index=which(DiagRips2$diagram[,1]==2)
    t2=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
    S[3,]=c(DiagRips2$diagram[t2,2],DiagRips2$diagram[t2,3])
  }else{
    t2=-1
    S[3,]=c(0,0)
  }
  
  # Combine the data.
  if(t1!=-1){
    data1=matrix(as.vector(DiagRips1[["cycleLocation"]][[t1]]),
                 ncol=3,byrow = F)
  }else{
    data1=c()
  }
  if(t2!=-1){
    data2=matrix(as.vector(DiagRips2[["cycleLocation"]][[t2]]),ncol=3,byrow = F)
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
  data=unique(data) # This step is to prevent the duplicated case and speed it up.
  
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
  
  death1=DiagRips1$deathLocation[t1,]
  death2=DiagRips2$deathLocation[t2,]
  data1=data1[data1[,1]!=-gap,]
  data2=data2[data2[,1]!=-gap,]
  #BirthRecal=max(BirthRecal(data1,data2),S[2,1],S[3,1]) # This birth estimate is for general.
  BirthRecal=max(S[2,1],S[3,1])
  seed=c(1,2,3,4) # Can be cahnged to random numbers.
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
