# This code is used to test how the combination behaving.
# Return as a matrix with full data birth time and death time and two
# cutting pieces' birth time and death times.
# n is the number of supplemental points.

# CUtting into half.

DiagCir <- function(X,gap,maxscale,maxdimension,n){
  X1=X[-X[,1]>gap,]
  X2=X[-X[,1]<gap,]
  Y=matrix(c(rep(-gap,n),seq(-2,2,4/(n-1))),nrow=n,byrow = F)
  S=matrix(0,nrow = 4,ncol = 2)
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
    X = rbind(X1,Y), maxdimension = maxdimension, maxscale = maxscale,
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
    X = rbind(X2,Y), maxdimension = maxdimension, maxscale = maxscale,
    library = "Dionysus", location = TRUE, printProgress = F)
  if( any(DiagRips2$diagram[,1]==1) ){
    index=which(DiagRips2$diagram[,1]==1)
    t2=index[which.max(DiagRips2$diagram[index,3]-DiagRips2$diagram[index,2])]
    S[3,]=c(DiagRips2$diagram[t2,2],DiagRips2$diagram[t2,3])
  }else{
    t2=-1
    S[3,]=c(0,0)
  }
  
  # Combine the data.
  if(t1!=-1){
    data1=matrix(as.vector(DiagRips1[["cycleLocation"]][[t1]]),ncol=2,byrow = F)
  }else{
    data1=c()
  }
  if(t2!=-1){
    data2=matrix(as.vector(DiagRips2[["cycleLocation"]][[t2]]),ncol=2,byrow = F)
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
  if( any(DiagRips3$diagram[,1]==1) ){
    index=which(DiagRips3$diagram[,1]==1)
    t3=index[which.max(DiagRips3$diagram[index,3]-DiagRips3$diagram[index,2])]
    S[4,]=c(DiagRips3$diagram[t3,2],DiagRips3$diagram[t3,3])
  }else{
    S[4,]=c(0,0)
  }
  return(S)
}