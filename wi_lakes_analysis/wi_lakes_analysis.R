# This script carries out the Wisconsin Lakes analysis in the manuscript by  
##Chenghui Li and Jessi Cisewski-Kehe titled, 
##"A Divide-and-Conquer Approach to Persistent Homology." 
#Full citation:  [To be added once available]

# Data source:  https://apps.dnr.wi.gov/lakes/lakepages/Results.aspx
# Cleaned data file used in paper:  lakes_wi_updated.csv
# DaC persistence diagrams:  

# Required packages
require(tidyverse)
require(this.path) #to autoselect working directory (ignore if you set manually)

# Working directory
setwd(this.path::here())
getwd() #Ensure path is to wi_lakes_analysis folder

# Load cleaned WI Lake data
lakes = read_csv("lakes_wi_updated.csv")

# Plot of lake locations
ggplot(lakes) +
  geom_point(aes(x=long,y=lat))

# Load WI lakes persistence diagrams
