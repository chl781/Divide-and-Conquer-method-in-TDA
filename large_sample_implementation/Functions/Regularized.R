# Regularize the list in Diag to prevent the trivial case appear.

Regularized <- function(S){
  if(is.vector(S$diagram)==T){
    names1=names(S$diagram)
    S$diagram=matrix(S$diagram,nrow=1)
    colnames(S$diagram)<-names1
    S$birthLocation=matrix(S$birthLocation,nrow=1)
    S$deathLocation=matrix(S$deathLocation,nrow=1)
    S$Doubt=matrix(S$Doubt,nrow=1)
  }
  return(S)
}