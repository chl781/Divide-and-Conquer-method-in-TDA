# My own calculation process for death time Eddited.
# x is the 
# death1 is the death time for x. death2 is the death time for the other data set.
# deathtime is the death time for x.

DeathRecal1<- function(data,death1,death2){
  a=0
  n=nrow(data)
  deathdistance=sqrt(sum((death1-death2)^2))# Add an extra estimates
  for(i in 1:n){
    for(j in 1:n){
      d1=MinLength(data[i,],data[j,],death1)
      d2=MinLength(data[i,],data[j,],death2)
      d=max(d1,d2)
      if(d>a){
        a=d
      }
    }
  }
  return(a) # death distance cannot be guaranteed to behave well.
}