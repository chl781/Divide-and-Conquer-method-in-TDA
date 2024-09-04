# 1 generate data for 100 times


j1=400

maxscale=2
maxdimension=1 # maxdimension = 2

for(i in 1:100){
  # Uniform generate data, and the distance between two circles is far.
  X1=rnorm(j1)
  X2=rnorm(j1)
  X <- cbind(X1,X2)
  X <- X/2/sqrt(rowSums(X^2))
  X<- X - matrix(c(-.1,.3),ncol=2,nrow=j1,byrow=T)
  
  X1=rnorm(j1/30)# 15
  X2=rnorm(j1/30)# 15
  Y <- cbind(X1,X2)
  Y <- Y/30/sqrt(rowSums(Y^2))
  Y <- Y+matrix(c(20/30,10/30),ncol=2,nrow=j1/30,byrow=T) #j1/15
  
  
  X<-rbind(X,Y)
  
  write.table( X,paste0("2closeCircles_",i,".csv"), sep=",",  col.names=FALSE,row.names = F)
  #write.csv(X,paste0("2closeCircles_",i,".csv"),row.names = F)
  
}





# 1 create figures for other datas

for (i in 1:100) {
  Dataname=paste0(i,"_output.csv")
  X=read.csv(Dataname, header = F)
  i=which(X[,1]==1)
  b_kmeans=X[i,2]
  d_kmeans=X[i,3]
  
}

