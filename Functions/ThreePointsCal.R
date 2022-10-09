# This is Jessi's guess. Draw 3 points out and calculate the max, min or average.

ThreePointsCal<- function(x,y){
  x1=NROW(x)
  y1=NROW(y)
  if(x1<=1){
    yDraw=y[sample(1:y1,3-x1),]
    return(rbind(x[1,],yDraw))
  } else if(y1<=1){
    xDraw=x[sample(1:x1,3-y1),]
    return(rbind(y[1,],xDraw))
  }else{
    xDraw=x[sample(1:x1,2),]
    yDraw=y[sample(1:y1,1),]
    return(rbind(xDraw,yDraw))
  }
}