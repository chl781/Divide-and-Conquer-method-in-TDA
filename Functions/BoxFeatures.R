# This is for ploting box based on the specific features.

# This function will plot the box based on the position of features.

BoxFeatures <- function(X,Diagcombine,...){
  plot(X,...)
  j=2
  for (i in 1:nrow(Diagcombine$diagram)) {
    lim1=Diagcombine$lim[[i]]
    box1=matrix(0,nrow=5,ncol=2)
    box1[1,]=lim1[,1]
    box1[2,]=c(lim1[1,1],lim1[2,2])
    box1[3,]=c(lim1[1,2],lim1[2,2])
    box1[4,]=c(lim1[1,2],lim1[2,1])
    box1[5,]=c(lim1[1,1],lim1[2,1])
    lines(box1,col=j)
    j=j+1
  }
}