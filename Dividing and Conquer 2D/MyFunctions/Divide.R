
# Generate the dividing diagrams' form.


# For now, assume it is for 2d problem.
Divide <- function(X,lim,n,by,m0){
  # Innitial S
  if(nrow(lim)!=ncol(X)){
    stop('The dimension for X or lim is wrong')
  }
  S=list()
  if(ncol(X)==2){
    Grid0 <- seq(lim[1,1], lim[1,2], by = by)
    Grid1 <- seq(lim[2,1], lim[2,2], by = by)
    Grid <- expand.grid(Grid0, Grid1)
    DTM <- dtm(X = X, Grid = Grid, m0 = m0)
    Xseq <- Grid0
    Yseq <- Grid0
    
    Xdtm=matrix(DTM, ncol = length(Yseq), nrow = length(Xseq))
    colnames(Xdtm)<-Yseq
    rownames(Xdtm)<-Xseq
    # n is the number it will break for each side.
    k=1
    for(i in 1:n){
      for(j in 1:n){
        # n is the number it will break for each side.
        stepsize1=(lim[1,2]-lim[1,1])/(n+1)
        stepsize2=(lim[2,2]-lim[2,1])/(n+1)
        lim1=c(lim[1]+(i-1)*stepsize1,lim[1]+(i+1)*stepsize1)
        lim2=c(lim[1]+(j-1)*stepsize2,lim[1]+(j+1)*stepsize2)
        Xdtmloc <- Xdtm[as.numeric(rownames(Xdtm))<=lim1[2]&as.numeric(rownames(Xdtm))>=lim1[1],
                        as.numeric(colnames(Xdtm))<=lim2[2]&as.numeric(colnames(Xdtm))>=lim2[1]]
        Diag1<-gridDiag(X,FUNvalues = Xdtmloc,lim=
                        matrix(c(lim1[1],lim1[2],lim2[1],lim2[2]),2),
                        by=by,m0=m0, sublevel = TRUE,
                        printProgress = TRUE,location=T,maxdimension = 2,library="Dionysus")
        Diag1$lim=rbind(lim1,lim2)
        S[[k]]=Diag1 # Add the diagram information into the list
        k=k+1
      }
    }
  }
  if(ncol(X)==3){
    
  }
  return(S)
}

