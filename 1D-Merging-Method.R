# This file loads a 2D data and generates 1d feature based on DaC method.

# Load packages

library(pacman)
p_load("dplyr","plotrix","spatstat","TDA","hitandrun","functional","Rfast","plotly","viridis","plot3D","ggplot2")
require(vrmlgen)
require(rgl)
require(fields)
require(knitr)
library("profmem") # used for testing memory
require(circular)
require(proxy)
require(wordspace)
require(geometry)


# Load function files
source('Functions/DiagCirSimp.R')# Add another Esimate method.
source('Functions/DiagCir4Pieces.R')
source('Functions/PlotRepeat.R')
source('Functions/PlotRepeat1.R')
source('Functions/BirthRecal2.R')
source('Functions/DeathRecal0.R')
source('Functions/DeathRecal1.R')
source('Functions/DeathRecal2.R')
source('Functions/DeathRecal_Circle.R')
source('Functions/ThreePointsCal.R')
source('Functions/MinLength.R')
source('Functions/DiagCir3d.R')
source('Functions/DeathRecal_Sphere.R')
source('Functions/Continuous.R')
source('Functions/DiagCirCont4.R')
source('Functions/DiagCont3d.R')
source('Functions/DiagContm2.R')
source('Functions/BoundaryConnect.R')
source('Functions/Utils.R')
source('Functions/Utils2.R')
source('Functions/Projected_Merge.R')
source('Functions/Diagm3Combine.R')
source('Functions/boundFind.R')
source('Functions/Matching2.R')

# Load 3D function files
source("Functions3Combine/3DiagContm2.R")
source("Functions3Combine/BoundaryConnect_.R")
source("Functions3Combine/Diagm3Combine_.R")
source("Functions3Combine/Utils_.R")
source("Functions3Combine/Utils2_.R")
source("Functions3Combine/BirthRecal2_.R")
source("Functions3Combine/DeathRecal0_.R")
source("Functions3Combine/DeathRecal1_.R")

# Parameter Setup
m=4 # Split the data into 7*7 subregions

range=matrix(c( -1, 1, -1, 1),nrow = 2, byrow=T)
maxscale=6

# Equal spaced sub-regions
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)

# Maxdimension setup
maxdimension=1

# Load data and do basic transform
j2=1
X=read.csv2(paste0("data/1Dcircle/2closeCircles_",j2,".csv"),header=F,sep=",")
X=as.numeric(as.matrix(X))
X=matrix(X,ncol=2)


# Generate divide data in sub-regions.
m=15
X_split=matrix(list(),m-1,m-1)
range<-matrix(c(-1-1.7/m,-1-1.4/m,.7+1.7/m,.4+1.4/m),2)
Matching_error=0.1

gap1=seq(range[1,1],range[1,2],length.out =m)
gap2=seq(range[2,1],range[2,2],length.out =m)

for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
    X_split[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
  }
}


maxdimension=1
maxscale=4
error=1

# Have the split diagrams
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range)
Combine1 = Diagm3Combine(X_split,m,Diag_split,range,maxdimension,maxscale,error)
PD = Combine1$diagram


# Have the split diagrams
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range)
Combine1 = Diagm3Combine(X_split,m,Diag_split,range,maxdimension,maxscale,error)
PD = Combine1$diagram

# This function works for d=1 case.

# Combine1 is the merged features.

# Plot for Combine1 representative data points.
plot3D::scatter2D(Combine1[[2]][[1]][,1],Combine1[[2]][[1]][,2])

# Diagram
Combine1$diagram
