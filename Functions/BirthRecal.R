# This program is to calculate the birth time of the union of a and b

BirthRecal <-function(a,b){
  s=c()
  t=0
  if(NROW(a[a[,2]>0,])>0&NROW(b[b[,2]>0,])>0){
    a1=which.max(a[a[,2]>0,1])
    b1=which.min(b[b[,2]>0,1])
    s=c(s,sqrt(sum((a[which(a[,2]>0)[a1],]-b[which(b[,2]>0)[b1],])^2)))
    t=t+1
  }
  if(NROW(a[a[,2]<0,])>0&NROW(b[b[,2]<0,])>0){
    a2=which.max(a[a[,2]<0,1])
    b2=which.min(b[b[,2]<0,1])
    s=c(s,sqrt(sum((a[which(a[,2]<0)[a2],]-b[which(b[,2]<0)[b2],])^2)))
    t=t+1
  }
  if(t<=1){
    warning('There is at most one side has data points.')
    return(0)
  }
  return(max(s))
}