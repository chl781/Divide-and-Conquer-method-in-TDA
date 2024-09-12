# My own calculation process for death time.
# x is the 
# death1 is the death time for x. death2 is the death time for the other data set.
# deathtime is the death time for x.

DeathRecal<- function(x,death1,death2, deathtime){
  for (i in nrow(x)) {
    if(sum((x[i,]-death1)^2)==deathtime^2)
      break
  }
  return(sqrt(sum((death2-x[i,])^2)))
}