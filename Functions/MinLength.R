

MinLength<-function(a,b,c){
  a1=sqrt(sum((a-b)^2))
  a2=sqrt(sum((c-b)^2))
  a3=sqrt(sum((a-c)^2))
  return(min(a1,a2,a3))
}