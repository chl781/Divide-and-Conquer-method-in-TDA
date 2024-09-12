
BirthRecal2_ <- function(Diag_split,Combined_diag_indices,ind_suspicious){
  estimate=0
  for(i in 1:NROW(Combined_diag_indices) ){
    idx1=Combined_diag_indices[i,]
    a1=as.numeric( Diag_split[[ idx1[1],idx1[2],idx1[3] ]]$diagram[
      ind_suspicious[[ idx1[1],idx1[2],idx1[3] ]][ idx1[4] ],2 ] )
    if(estimate<a1 ){
      estimate=a1
    }
  }
  
  return(estimate)
}