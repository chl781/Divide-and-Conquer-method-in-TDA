# This script considers a 2D point cloud (a big circle and a little circle)
## and generates an H1 persistence diagram based on DaC method.



# Load packages
library(tidyverse)
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
source('Functions/DiagCirSimp.R')# Add another Estimate method.
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

# Input parameters for DaC and persistent homology
maxdimension=1 # Maximum homology dimension considered
maxscale=4 # Maximum distance scale considered
error=1  # [Chenghui: add brief explanation for this]

# Load data and do basic transform 
which_data_set = 5 # Select which iid data set from {1, 2, ..., 100}
X=read.csv(paste0("data/1Dcircles/2closeCircles_",which_data_set,".csv"),header=FALSE)
X = as.matrix(X, ncol = 2)

# Quick visualization of data (one big circle and one small circle)
ggplot(as_tibble(X), aes(x=V1, y=V2)) +
  geom_point(size=5) +
  labs(x="X",y="Y") +
  theme_bw()

# Generate divide data in sub-regions.
m=15 # Number of subregions
X_split=matrix(list(),m-1,m-1)
range<-matrix(c(-1-1.7/m,-1-1.4/m,.8+1.4/m,.4+1.4/m),2) # [Chenghui: add brief explanation for this]
Matching_error=0.1 # [Chenghui: add brief explanation for this]

gap1=seq(range[1,1],range[1,2],length.out =m)
gap2=seq(range[2,1],range[2,2],length.out =m)


# [Chenghui: add brief explanation for this]
for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
    X_split[[i,j]]=X[X[,1]>gap1[i]&X[,2]>gap2[j]&X[,1]<gap1[i+1]&X[,2]<gap2[j+1],]
  }
}


# Have the split diagrams  # [Chenghui: what does this mean??]
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range) # [Chenghui: what is this doing?]
Combine1 = Diagm3Combine(X_split,m,Diag_split,range,maxdimension,maxscale,error) # [Chenghui: what is this doing?] # Combine1 is the merged features.
PD = Combine1$diagram # [Chenghui: what is this?  Final DaC persistence diagram?]
# This function works for d=1 case.



# Plot for Combine1 representative data points.
plot3D::scatter2D(Combine1[[2]][[1]][,1],Combine1[[2]][[1]][,2]) # [Chenghui: what is this?]
plot3D::scatter2D(Combine1[[2]][[2]][,1],Combine1[[2]][[2]][,2]) # [Chenghui: what is this?]

# H1 DaC Diagram
dac_diagram_h1 = Combine1$diagram

# True H1 Diagram based on full data using TDA package
true_diagram = ripsDiag(X, maxdimension = maxdimension, maxscale = maxscale)
true_diagram_h1 = true_diagram$diagram[true_diagram$diagram[,1]==1, ] #True H1 Diagram

# Compare DaC and true H1 diagrams
true_diagram_h1[order(true_diagram_h1[, 2]),]
dac_diagram_h1[order(dac_diagram_h1[, 2]),]

