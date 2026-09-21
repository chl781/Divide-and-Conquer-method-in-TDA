# Small, median and large features recovery

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

# Read data from write.csv(circle_data,file = "disjoint_circle_manifolds.csv",row.names = FALSE)
circle_data <- read.csv("disjoint_circle_manifolds.csv",row.names = FALSE)


# DaC method

##Load data file
circle_data = read.csv("data/Disjoint_Manifolds/disjoint_circle_manifolds.csv")
  
  
# Parameter Setup

m=8
X=circle_data[,1:2]

range=matrix(c( min(X[,1]), max(X[,1]),min(X[,2]), max(X[,2])),nrow = 2, byrow=T)
maxscale=5

# Equal spaced sub-regions
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)

# Maxdimension setup
maxdimension=1

# Load data and do basic transform
error=0.1

# Generate divide data in sub-regions.
X_split = array(list(),c(m-1,m-1))
for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
      X_split[[i,j]] = X[ X[,1] >= gap1[i] & X[,1] < gap1[i+1] 
                            & X[,2] >= gap2[j] & X[,2] < gap2[j+1],]
  }
}




# Have the split diagram - this may take a couple minutes to run
Diag_split <- DiagContm2(X,m,maxscale,maxdimension,range)
Combine1 = Diagm3Combine(X_split,m,Diag_split,range,maxdimension,maxscale,error)
PD = Combine1$diagram



# Plot the diagram
birth <- PD[, 2]
death <- PD[, 3]

lims <- range(c(birth, death), finite = TRUE)


df <- data.frame(birth = birth, death = death)

ggplot(df, aes(x = birth, y = death)) +
  geom_point(
    size  = 6,
    color = "black",
    shape = 16
  ) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(1) +
  xlab("Birth") +
  ylab("Death") +
  xlim(lims) +
  ylim(lims) +
  theme_bw() +
  theme(
    text = element_text(size = 40),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )
