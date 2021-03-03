
# Determine if there is features touch the boundary
# tol is the level we can treat as close to the boundary
# For now, we set tol as 0 to test.

IdentifyBoun <- function(D,lim,tol){
  Result=1
  if(any(abs(D[,1,1]-lim[1,1])<=tol|abs(D[,1,1]-lim[1,2])<=tol|
         abs(D[,1,2]-lim[2,1])<=tol|abs(D[,1,2]-lim[2,2])<=tol)){
    Result=0
  }
  return(Result)
}