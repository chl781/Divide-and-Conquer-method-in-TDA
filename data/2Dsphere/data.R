# Generate some spheres data
set.seed(2023)

j1=1600
n=j1

for(i in 1:100){
  X1=rnorm(n)
  X2=rnorm(n)
  X3=rnorm(n)
  
  X=cbind(X1,X2,X3)
  X=X/sqrt(rowSums(X^2))
  write.csv2(X,paste0("data",i,".csv"),col.names = F,row.names=F)
}