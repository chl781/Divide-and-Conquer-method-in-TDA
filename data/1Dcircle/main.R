
# Load packages

library(pacman)
p_load("dplyr","plotrix","spatstat","TDA","hitandrun","functional","Rfast","plotly","viridis","plot3D","ggplot2")
require(rgl)
require(fields)
require(knitr)
library("profmem") # used for testing memory
require(circular)
require(proxy)
require(wordspace)
require(geometry)


# Load function files
source('../../Functions/DiagCirSimp.R')
source('../../Functions/DiagCir4Pieces.R')
source('../../Functions/PlotRepeat.R')
source('../../Functions/PlotRepeat1.R')
source('../../Functions/BirthRecal2.R')
source('../../Functions/DeathRecal0.R')
source('../../Functions/DeathRecal1.R')
source('../../Functions/DeathRecal2.R')
source('../../Functions/DeathRecal_Circle.R')
source('../../Functions/ThreePointsCal.R')
source('../../Functions/MinLength.R')
source('../../Functions/DiagCir3d.R')
source('../../Functions/DeathRecal_Sphere.R')
source('../../Functions/Continuous.R')
source('../../Functions/DiagCirCont4.R')
source('../../Functions/DiagCont3d.R')
source('../../Functions/DiagContm2.R')
source('../../Functions/BoundaryConnect.R')
source('../../Functions/Utils.R')
source('../../Functions/Utils2.R')
source('../../Functions/Projected_Merge.R')
source('../../Functions/Diagm3Combine.R')
source('../../Functions/boundFind.R')
source('../../Functions/Matching2.R')

# Load 3D function files
source("../../Functions3Combine/3DiagContm2.R")
source("../../Functions3Combine/BoundaryConnect_.R")
source("../../Functions3Combine/Diagm3Combine_.R")
source("../../Functions3Combine/Utils_.R")
source("../../Functions3Combine/Utils2_.R")
source("../../Functions3Combine/BirthRecal2_.R")
source("../../Functions3Combine/DeathRecal0_.R")
source("../../Functions3Combine/DeathRecal1_.R")

# Parameter Setup
m=5 # Split the data into 4*4*4 subregions

range=matrix(c( -1, 1,-1, 1),nrow = 2, byrow=T)
maxscale=6

# Equally spaced sub-regions
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)

# Maxdimension setup
maxdimension=1

# Load data and do basic transform
j2=1
X=read.csv2(paste0("data",j2,".csv"),header=T,sep=";")
X=as.numeric(as.matrix(X))
X=matrix(X,ncol=2)


# Generate divide data in sub-regions.
X_split = array(list(),c(m-1,m-1))
for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
    X_split[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
  }
}

error=0.1 # Define the matching error bound for the cancellation method.


# Have the split diagram
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range)
Combine1 = Diagm3Combine(X_split,m,Diag_split,range,maxdimension,maxscale,error)
PD = Combine1$diagram

