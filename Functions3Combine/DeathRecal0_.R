DeathRecal0_<- function(data){
  a=0
  n=NROW(data)
  c=0
  DeathRecal1=0
  #alocation=matrix(0,nrow=2,ncol=3)
  for(i in 1:(n-3)){
    for(j in (i+1):(n-2)){
      a1=sqrt(sum((data[i,]-data[j,])^2))
      for(k in (j+1):(n-1)){
        a2=sqrt(sum((data[i,]-data[k,])^2))
        a3=sqrt(sum((data[j,]-data[k,])^2))
        for (l in (k+1):n) {
          a4=sqrt(sum((data[i,]-data[l,])^2))
          a5=sqrt(sum((data[j,]-data[l,])^2))
          a6=sqrt(sum((data[k,]-data[l,])^2))
          a= min(c(a1,a2,a3,a4,a5,a6))
          if(a>DeathRecal1){
            DeathRecal1=a
            #alocation=rbind(data[i,],data[k,],data[j,])
            c= max(c(a1,a2,a3,a4,a5,a6))
          }
        }
      }
    }
  }
  return(c) # death distance can be guaranteed to approximate well.
}

