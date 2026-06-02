# Stanford Bunny dataset
# This file is not for the csv but for the original 3D data preprecessing and running DaC method.

#install.packages("Rvcg")
#install.packages("rgl")

library(Rvcg)
library(rgl)

# Load the bunny mesh from original file.
# The file can be downloaded from the Stanford 3D Scanning Repository: http://graphics.stanford.edu/data/3Dscanrep/
mesh <- vcgPlyRead("bun_zipper.ply", updateNormals = TRUE)

# mesh$vb contains vertices in homogeneous coordinates (4 × n matrix)
points <- t(mesh$vb[1:3, ])   # Convert to n×3 matrix

head(points)

# Visualize the dataset
plot3d(points, size = 3, col = "gray")

### The following is a preprocessing step.

mesh_ds <- vcgUniformRemesh(
  mesh,
  voxelSize = 0.004,   
  discretize = TRUE,
  multiSample = FALSE
)

pc_down <- t(mesh_ds$vb[1:3, ])

nrow(pc_down)
plot3d(pc_down, size = 4, col = "gray")


# Save pc_down to csv
output_filename = "Stanford_Bunny.csv"
write.table(
  pc_down,
  file = output_filename,
  sep = ",",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# Load packages for DaC method and visualization.
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
# Split the data into 7*7*7 subregions
m=8

range=matrix(c( min(pc_down[,1]), max(pc_down[,1]),min(pc_down[,2]), max(pc_down[,2]), min(pc_down[,3]), max(pc_down[,3])),nrow = 3, byrow=T)
maxscale=1

# Equal spaced sub-regions
gap1=seq(range[1,1], range[1,2], length.out = m)
gap2=seq(range[2,1], range[2,2], length.out = m)
gap3=seq(range[3,1], range[3,2], length.out = m)

# Maxdimension setup
maxdimension=2

# Load data and do basic transform
j2=1
X=pc_down


# Generate divide data in sub-regions.
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




# Have the split diagram
Diag_split=DiagContm2_(X,m,X_split,gap1,gap2,gap3,maxscale,maxdimension) # This function works for d=2 case.

# Merge the sub-regions
Combine1=Diagm3Combine_(X_split,m,Diag_split,
                        gap1,gap2,gap3,maxdimension,maxscale)

# Combine1 is the merged features.

# Plot for Combine1 representative data points.
plot3D::scatter3D(Combine1[[2]][[1]][,1],Combine1[[2]][[1]][,2],Combine1[[2]][[1]][,3])

# Diagram
Combine1$diagram


## Visualize the result

library(ggplot2)

df <- data.frame(
  birth = Combine1$diagram[,"birth"],
  death = Combine1$diagram[,"death"]
)


# expand limits (increase x-axis scale safely)
lims <- c(0, max(df$death) * 1.25)


ggplot(df, aes(birth, death)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = lims, expand = c(0, 0)) +
  scale_y_continuous(limits = lims, expand = c(0, 0)) +
  coord_fixed(1) +
  labs(
    x = "Birth",
    y = "Death",
    title = expression("Persistence Diagram (" * H[2] * ")")
  ) +
  theme_minimal(base_size = 14) +   # global base size
  theme(
    plot.title = element_text(size = 25, face = "bold"),
    axis.title = element_text(size = 30, face = "bold"),
    axis.text  = element_text(size = 20)
  )



# Save the data illustration 

library(rgl)

# Open off-screen device (important for saving)
open3d(useNULL = FALSE)

# Plot the point cloud
plot3d(pc_down,
       size = 4,
       col = "gray30",
       decor = FALSE)   # <-- removes axes, ticks, labels, box

# Clean background
bg3d(color = "white")

# Lighting
light3d(theta = 30, phi = 30, viewpoint.rel = TRUE)
light3d(theta = -30, phi = -30, viewpoint.rel = TRUE)

# View
view3d(theta = 40, phi = 20, zoom = 0.75)

# High resolution window
par3d(windowRect = c(100, 100, 1400, 1400))

# Save
rgl.snapshot("bunny_data.png", top = TRUE)


## In the following, we report the memory cost
library(peakRAM)

peakRAM(
  Diag_split <- DiagContm2_(X, m, X_split,
                            gap1, gap2, gap3,
                            maxscale, maxdimension)
)

peakRAM(
  Combine1<-Diagm3Combine_(X_split,m,Diag_split,
                                gap1,gap2,gap3,maxdimension,maxscale))
