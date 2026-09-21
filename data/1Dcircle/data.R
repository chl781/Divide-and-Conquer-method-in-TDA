# Generate some spheres datadata from 1D circle.
set.seed(2023)

j1=400
n=j1

for(i in 1:100){
  X1=rnorm(n)
  X2=rnorm(n)
  
  X=cbind(X1,X2)
  X=X/sqrt(rowSums(X^2))
  write.csv2(X,paste0("data",i,".csv"),col.names = F,row.names=F)
}