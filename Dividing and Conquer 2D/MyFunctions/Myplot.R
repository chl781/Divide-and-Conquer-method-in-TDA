# Plot including suspicious features.

Myplot <- function(x,...){
  if(is.vector(x$diagram)==T){
    names1=names(x$diagram)
    x$diagram=matrix(x$diagram,nrow=1)
    colnames(x$diagram)<-names1
  }
  if(x$Doubt[which.max(x$diagram[,'Death'])]==1){
    plot.diagram(x = x[["diagram"]][x$Doubt==1,],col=2,diagLim = c(0,max(x$diagram[,'Death'])),
                 ...)
    if(length(x$Doubt==0)>0){
      plot.diagram(x = x[["diagram"]][x$Doubt==0,],
                   col=3,add = T,diagLim = c(0,max(x$diagram[,'Death'])))
    }
  } else{
    plot.diagram(x = x[["diagram"]][x$Doubt==0,],col=3,diagLim = c(0,max(x$diagram[,'Death'])),
                 ...)
    if(length(x$Doubt==1)>0){
      plot.diagram(x = x[["diagram"]][x$Doubt==1,],
                   col=2,add = T,diagLim = c(0,max(x$diagram[,'Death'])))
    }
  }
}