# Create illustrating boxes for the dividing steps

HighLightBox <- function(X,Diag,...){
  plot(X,...)
  j=1
  for (i in Diag) {
    lim1=i$lim
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