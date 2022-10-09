# Combine different features by using the death location
# and then identify the indexes.

Death <- function(S){
  index=1
  A=S$diagram[1,]
  DeathLocation=S$deathLocation[1,]
  
  if(nrow(S$diagram)==1) return(index)
  for(i in 2:nrow(S$diagram)){
    if(is.vector(DeathLocation)==T){
      if(isTRUE(all.equal(DeathLocation,S$deathLocation[i,] ))==T){
        former=S$diagram[index,'Death']-S$diagram[index,'Birth']
        compete=S$diagram[i,'Death']-S$diagram[i,'Birth']
        if(former<compete){
          index=i
        }
      }else{
        DeathLocation=rbind(DeathLocation,S$deathLocation[i,])
        index=c(index,i)
      }
    }else{
      if(any(apply(DeathLocation,1,function(x, want) 
        isTRUE(all.equal(x, want)), S$deathLocation[i,]))){
        indexformer=which(apply(DeathLocation,1,function(x, want) 
                   isTRUE(all.equal(x, want)), S$deathLocation[i,])==1)
        former=S$diagram[indexformer,'Death']-S$diagram[indexformer,'Birth']
        compete=S$diagram[i,'Death']-S$diagram[i,'Birth']
        if(former<compete){
          index[indexformer]=i
        }
      } else{
        DeathLocation=rbind(DeathLocation,S$deathLocation[i,])
        index=c(index,i)
      }
    }
  }
  return(index)
}