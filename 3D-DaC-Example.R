# This script considers a 3D point cloud (sphere)
## and generates an H2 persistence diagram based on DaC method.

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


# Input parameters for DaC and persistent homology
maxdimension=2 # Maximum homology dimension considered
maxscale=6 # Maximum distance scale considered
m = 8 # Split the data into 7*7*7 subregions


# [Chenghui:what is this for?]
range=matrix(c( -1, 1,-1, 1, -1, 1),nrow = 3, byrow=T)


# The following defines equally spaced sub-regions
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)
gap3=seq(range[3,1], range[3,2], length.out = m)


# Load data and do basic transform
which_data_set = 5 # Select which iid data set from {1, 2, ..., 100}
X=read.csv2(paste0("data/2Dsphere/data",which_data_set,".csv"),header=T,sep=";")
X=as.numeric(as.matrix(X))
X=matrix(X,ncol=3)

# Quick visualization of data (sphere)
plot3D::scatter3D(X[,1],X[,2],X[,3], col=1)


# Partition data into subregions.
X_split = array(list(),c(m-1,m-1,m-1))
for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
    for(k in 1:(m-1)){
      X_split[[i,j,k]] = X[ X[,1] >= gap1[i] & X[,1] < gap1[i+1] 
                            & X[,2] >= gap2[j] & X[,2] < gap2[j+1]
                            &X[,3] >= gap3[k] & X[,3] < gap3[k+1],]
    }
  }
}



# Have the split diagram  # [Chenghui: what does this mean??] (This function works for H2 case with 3D data.)
Diag_split=DiagContm2_(X,m,X_split,gap1,gap2,gap3,maxscale,maxdimension) 

# Merge the subregions (This function works for H2 case with 3D data.)
Combine1=Diagm3Combine_(X_split,m,Diag_split, # Combine1 is the merged features.
                        gap1,gap2,gap3,maxdimension,maxscale)


# Plot for Combine1 representative data points.
plot3D::scatter3D(Combine1[[2]][[2]][,1],Combine1[[2]][[2]][,2],Combine1[[2]][[2]][,3]) 
# [Chenghui:  I changed this from Combine1[[2]][[1]] to Combine1[[2]][[2]], because the first feature had lower persistence]

# H2 DaC Diagram
dac_diagram_h2 = Combine1$diagram

# True H2 Diagram based on full data is (very) slow and may not run on a personal computer
##so the approximate "true" H2 diagram is based on a random subsample of the data X
## Note:  the following is slow, even for small n_sub
n_sub = 300 # Size of subsample
paste(n_sub/nrow(X)*100, "% of data used to compute approximate H2 diagram")
sub_X = X[sample(1:nrow(X),n_sub, replace=FALSE), ] # Sub-samble

library(ripserr) # This package computes diagrams for Vietoris-Rips filtration
approx_diagram = vietoris_rips(sub_X, max_dim = 2) # This is slow even for small sample due to max_dim = 2
## If you change max_dim to 1, it is much faster.
approx_diagram_h2 = cbind(approx_diagram$dimension, approx_diagram$birth, approx_diagram$death)[approx_diagram$dimension==2,]


# Compare DaC and (very approximate) true H2 diagrams
## Note:  These diagrams are not really comparable since they are based on different sample sizes
approx_diagram_h2[order(true_diagram_h2[, 2]),] # (Very) approximate H2 diagram based on small subsample
dac_diagram_h2[order(dac_diagram_h2[, 2]),] # DaC H2 diagram



