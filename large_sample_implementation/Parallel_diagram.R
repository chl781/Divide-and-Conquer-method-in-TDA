# Simulation_Only uniform distribution on Sphere. 
# That is the underlying manifold can only be perfect circles
# or perfect spheres. Then this problem can be reduced based on 1 method of calculating
# the death time.

# Consider generating some data from this kind of manifolds and do the analysis
# This is parallel case used in CHTC.

# Only generate the data following a wrappednormal distribution.

# So, the only change should be the distribution of the data. However, this should not 
# play a role.

# Number of data also play a role.

# In this file, only the data generation method is modified.

(args = commandArgs(trailingOnly=TRUE))
(j <- args[1])
j2 <- as.numeric(j)
j2 <- j2 + 1

j <- args[2]
j1 <- as.numeric(j)

###### Derive a MSE data.frame but for only one iteration.

# Not use the same random seed!

library(pacman)
p_load("dplyr","plotrix","spatstat","TDA","hitandrun","functional","Rfast","plotly","viridis","plot3D","ggplot2")
require(vrmlgen)
require(rgl)
require(fields)
require(knitr)
library("profmem") # used for testing memory
library(mph) # Linearized approximation to Rips Filtration
require(circular)

require(proxy)
require(wordspace)
require(geometry)


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

source("Functions3Combine/3DiagContm2.R")
source("Functions3Combine/BoundaryConnect_.R")
source("Functions3Combine/Diagm3Combine_.R")
source("Functions3Combine/Utils_.R")
source("Functions3Combine/Utils2_.R")

source("Functions3Combine/BirthRecal2_.R")
source("Functions3Combine/DeathRecal0_.R")
source("Functions3Combine/DeathRecal1_.R")

source("Functions/parallel/Utils2_.R")
source("Functions/parallel/DiagContm2_.R")


############################# Calculation.



X=read.csv2("data/lakes_wi_updated.csv",header=T,sep = ",")

df<- X[,c("lat","long")]
df_clean <- na.omit(df)

colnames(df_clean)=NULL
df_clean= matrix(as.numeric(as.matrix(df_clean)),ncol=2)


range=matrix(c( min(df_clean[,1]), max(df_clean[,1]),
                min(df_clean[,2]), max(df_clean[,2])),nrow = 2, byrow=T)

m=17
maxscale=max(df_clean[,1])+max(df_clean[,2])-min(df_clean[,1])-min(df_clean[,2])
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)

X_split = matrix(list(),m-1,m-1)

for (i in 1:(m-1)) {
  for(j in 1:(m-1)){
      X_split[[i,j]] = df_clean[ df_clean[,1] >= gap1[i] & df_clean[,1] < gap1[i+1] 
                            & df_clean[,2] >= gap2[j] & df_clean[,2] < gap2[j+1],]
      }
}


maxdimension=1
error=1

i=((j2-1) %/% (m-1)) + 1
j=(j2-1) %% (m-1) +1

X_process=X_split[[i,j]]

# Have the split diagrams
Diag_split=DiagContm2_(X_process,m,i,j,gap1,gap2)


folder_path="out/diag/"
saveRDS(Diag_split, paste0(folder_path,"diag_", i, "_",j,".rds"))



