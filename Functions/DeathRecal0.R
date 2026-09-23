# A substitute for DeathRecal1 to compute the death time for a data set. This is a more efficient method to compute the death time.

DeathRecal0<- function(data){
  a=0
  n=NROW(data)
  c=0
  DeathRecal1=0
  for(i in 1:(n-2)){
    for(j in (i+1):(n-1)){
      a1=sqrt(sum((data[i,]-data[j,])^2))
      for(k in (j+1):n){
        a2=sqrt(sum((data[k,]-data[j,])^2))
        a3=sqrt(sum((data[i,]-data[k,])^2))
        a= min(a1,a2,a3)
        if(a>DeathRecal1){
          DeathRecal1=a
          #alocation=rbind(data[i,],data[k,],data[j,])
          c= max(a1,a2,a3)
        }
      }
    }
  }
  return(c) # death distance can be guaranteed to approximate well.
}

