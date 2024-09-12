# This is another way to capture the Death time.

DeathRecal2<- function(data){
  a=0
  n=NROW(data)
  if(n==0) {return(0)}
  for(i in 1:n){
    for(j in 1:n){
      #if(min(sqrt(rowSums((data-
      #                    matrix(rep((data[i,]+data[j,])/2,n),
      #                                             ncol=2,byrow = T))^2)))==
      #   sqrt(sum((data[j,]-data[i,])^2))/2){
      #  if(sqrt(sum((data[j,]-data[i,])^2))>a)
      #  a=sqrt(sum((data[j,]-data[i,])^2))
      #}
      ################ Comment out the "if" sentence (Can be concave case?)

      b=  2*min(sqrt(rowSums((data-matrix(rep((data[i,]+data[j,])/2,n),
                              ncol=2,byrow = T))^2)))
      if(b>a)
        a=b
    }
  }
  return(a) # death distance cannot be guaranteed to behave well.
}