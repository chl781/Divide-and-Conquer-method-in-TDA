# Back up code

DeathRecal1_<- function(data,nsample=10^5){
  a=0
  n=NROW(data)
  c=0
  sampleseed=1:n
  DeathRecal1=0
  for(ij in 1:nsample){
	ijkl=sample(sampleseed,4,replace=F)
	i=ijkl[1]
	j=ijkl[2]
	k=ijkl[3]
	l=ijkl[4]
	a1=sqrt(sum((data[i,]-data[j,])^2))
        a2=sqrt(sum((data[i,]-data[k,])^2))
        a3=sqrt(sum((data[j,]-data[k,])^2))
        a4=sqrt(sum((data[i,]-data[l,])^2))
        a5=sqrt(sum((data[j,]-data[l,])^2))
        a6=sqrt(sum((data[k,]-data[l,])^2))
        a= min(c(a1,a2,a3,a4,a5,a6))
	if(a>DeathRecal1){
            	DeathRecal1=a
            	c= max(c(a1,a2,a3,a4,a5,a6))
          	}
	}


  return(c) # death distance can be guaranteed to approximate well.
}

