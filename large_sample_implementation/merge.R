(args = commandArgs(trailingOnly=TRUE))
(j <- args[1])
j2 <- as.numeric(j)
j2 <- j2 + 1

j <- args[2]
j1 <- as.numeric(j)

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
source("Functions3Combine/DiagContm2_.R")


############################# Calculation.

dist_bound=readRDS("data/dist_bound.rds")

X=read.csv2("data/lakes_wi_updated.csv",header=T,sep = ",")

df<- X[,c("lat","long")]
df_clean <- na.omit(df)

colnames(df_clean)=NULL
df_clean= matrix(as.numeric(as.matrix(df_clean)),ncol=2)


range=matrix(c( min(df_clean[,1]), max(df_clean[,1]),
                min(df_clean[,2]), max(df_clean[,2])),nrow = 2, byrow=T)

m=17
maxscale=max(df_clean[,1])+max(df_clean[,2])-min(df_clean[,1])-min(df_clean[,2])


# Have the split diagrams
Diag_suspicious <- ripsDiag(dist_bound, maxdimension=1, maxscale = maxscale,
                           library = "Dionysus", dist="arbitrary",
                           location = T)

folder_path="out/"
saveRDS(Diag_suspicious, paste0(folder_path,"diag_suspicious",".rds"))
