# An example with more circles.

## Generate the data.
n=500
X=rnorm(n)
Y=rnorm(n)

X=matrix(c(X,Y),nrow=n,ncol=2)
X1=X/sqrt(rowSums(X^2))

n=200
X=rnorm(n)
Y=rnorm(n)

X=matrix(c(X,Y),nrow=n,ncol=2)
X2=X/sqrt(rowSums(X^2))*1.5+c(3,3)

n=200
X=rnorm(n)
Y=rnorm(n)

X=matrix(c(X,Y),nrow=n,ncol=2)
X3=X/sqrt(rowSums(X^2))/2+c(-2,-2)

X=rbind(X1,X2,X3)


## Split the data and generate the sub-features.
m=8
X_split=matrix(list(),m-1,m-1)
#range<-matrix(c(-1,-1,1,1),2)
range <- matrix(c(-3,-3,5,5),2)
Matching_error=1

gap1=seq(range[1,1],range[1,2],length.out =m)
gap2=seq(range[2,1],range[2,2],length.out =m)

for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
    X_split[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
  }
}

maxdimension=1
maxscale=100

# Have the split diagrams
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range)

