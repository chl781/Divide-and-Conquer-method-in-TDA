# Combine Diagrams

# Input: Diagram list

CombineDiagram <- function(D,tol){
  # D should be a list
  if(is.list(list(1,2))==F){
    stop("Input needs to be diagram list.")
  }
  # Generate the innitial form.
  S=D[[1]]
  S$lim=rep(list(S$lim),nrow(S$diagram))
  # Combine all the information first
  for(i in D){
    if(i$lim[1,1]==D[[1]]$lim[1,1]&i$lim[1,2]==D[[1]]$lim[1,2]&
       i$lim[2,2]==D[[1]]$lim[2,2]&i$lim[2,1]==D[[1]]$lim[2,1]) next()
    S$diagram=rbind(S$diagram,i$diagram)
    S$birthLocation=rbind(S$birthLocation,i$birthLocation)
    S$deathLocation=rbind(S$deathLocation,i$deathLocation)
    S$cycleLocation=c(S$cycleLocation,i$cycleLocation)
    S$lim=c(S$lim,rep(list(i$lim),nrow(i$diagram)))
  }
  
  # Get rid of dimension 0 features.
  i=1
  if(NROW(S$diagram)==0){stop('no features detected!')}
  #NROW and NCOL is to prevent trivial result like nrow=1
  #ifelse(test=(max(NROW(S$diagram),NCOL(S$diagram)) >= 3), yes=max(NROW(S$diagram),NCOL(S$diagram)), no=min(NROW(S$diagram),NCOL(S$diagram)))
  while(i<=ifelse(test=(max(NROW(S$diagram),NCOL(S$diagram)) > 3), yes=max(NROW(S$diagram),NCOL(S$diagram)), no=min(NROW(S$diagram),NCOL(S$diagram)))){
    if(is.vector(S$diagram)==T){
      names1=names(S$diagram)
      S$diagram=matrix(S$diagram,nrow=1)
      colnames(S$diagram)=names1
    }else if(as.matrix(S$diagram)[i,'dimension']==0){
      S$diagram=S$diagram[-i,]
      S$birthLocation=S$birthLocation[-i,]
      S$deathLocation=S$deathLocation[-i,]
      S$cycleLocation[[i]]=NULL
      S$lim[[i]]=NULL
      i=i-1
    }
    i=i+1
  }
  if(NROW(S$diagram)==0){stop('no features detected!')}
  S$Doubt=rep(1,ifelse(test=(max(NROW(S$diagram),NCOL(S$diagram)) > 3), yes=max(NROW(S$diagram),NCOL(S$diagram)), no=min(NROW(S$diagram),NCOL(S$diagram))))
  # Innitialize doubt labels.
  # If doubt is 1, then it is undoubtful, otherwise 0 it is doubtful.
  for(j in 1:ifelse(test=(max(NROW(S$diagram),NCOL(S$diagram)) > 3), yes=max(NROW(S$diagram),NCOL(S$diagram)), no=min(NROW(S$diagram),NCOL(S$diagram)))){
    S$Doubt[j]=IdentifyBoun(S$cycleLocation[[j]],S$lim[[j]],tol=tol)
  }
  
  # We will need to test if there is an identity for the death location.
  # We only leave the ones having the largest significance level
  # and use it to identify if it is suspicious
  
  S=Regularized(S)
  
  s=Death(S)
  
  
  # Save the important features.
  if(length(s)>=1){
    S$diagram=S$diagram[s,]
    S$birthLocation=S$birthLocation[s,]
    S$deathLocation=S$deathLocation[s,]
    S$Doubt=S$Doubt[s]
    S$cycleLocation=S$cycleLocation[s]
    S$lim=S$lim[s]
  }
  
  S=Regularized(S)
  return(S)
}
