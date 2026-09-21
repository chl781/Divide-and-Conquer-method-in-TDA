# 1 generate data for 100 times


j1=400

maxscale=2
maxdimension=1 

for(i in 1:100){
  # Uniform generate data, and the distance between two circles is far.
  X1=rnorm(j1)
  X2=rnorm(j1)
  X <- cbind(X1,X2)
  X <- X/2/sqrt(rowSums(X^2))
  X<- X - matrix(c(-.1,.3),ncol=2,nrow=j1,byrow=T)
  
  X1=rnorm(j1/30)
  X2=rnorm(j1/30)
  Y <- cbind(X1,X2)
  Y <- Y/30/sqrt(rowSums(Y^2))
  Y <- Y+matrix(c(20/30,10/30),ncol=2,nrow=j1/30,byrow=T)
  
  
  X<-rbind(X,Y)
  
  write.table( X,paste0("2closeCircles_",i,".csv"), sep=",",  col.names=FALSE,row.names = F)
  #write.csv(X,paste0("2closeCircles_",i,".csv"),row.names = F)
  
}



