# This script carries out the Wisconsin Lakes analysis in the manuscript by  
##Chenghui Li and Jessi Cisewski-Kehe titled, 
##"A Divide-and-Conquer Approach to Persistent Homology." 
#Full citation:  [To be added once available]

# Data source:  https://apps.dnr.wi.gov/lakes/lakepages/Results.aspx
# Cleaned data file used in paper:  lakes_wi_updated.csv
# DaC persistence diagrams:  

# Required packages
require(tidyverse)
require(TDA)
require(this.path) #to autoselect working directory (ignore if you set manually)

# Working directory
setwd(this.path::here())
getwd() #Ensure path is to wi_lakes_analysis folder

# Load persistence diagram plotting function
source("GetDiagram.R")


################################### Load and visualize lake data
# Load cleaned WI Lake data
lakes = read_csv("lakes_wi_updated.csv")
partition_coordinates = readRDS("partition_coordinates.rds") #lat/long for southern and northern partitions

# Plot of lake locations
ggplot(lakes) +
  geom_point(aes(x=long,y=lat))


################################### Load and visualize persistence diagrams
# Load WI lakes persistence diagrams (H1 features only!)
pd_north = readRDS("pd_north.rds") #diagrams for northern region 
pd_south = readRDS("pd_south.rds") #diagrams for southern region 


# Visualize one of the persistence diagrams 
i=5 #Select i in {1,2,...,8}
GetDiagram(pd_north[[i]],labels=c(expression(H[1])))
GetDiagram(pd_south[[i]],labels=c(expression(H[1])))



################################### Compute landscape functions
num_layers=4 # Set number of landscape function layers to retain 
land_north <- rep(list(matrix(NA, ncol=8, nrow=1000)),num_layers) #Store northern landscapes
land_south <- rep(list(matrix(NA, ncol=8, nrow=1000)),num_layers) #Store southern landscapes

# Define landscape function domain (tseq) based on H1 persistences:
tmin = min(sapply(1:8, function(ii) min(pd_north[[ii]][,2:3], pd_south[[ii]][,2:3])))
tmax = max(sapply(1:8, function(ii) max(pd_north[[ii]][,2:3], pd_south[[ii]][,2:3])))
tseq = seq(tmin, tmax, length=1000)

for(land_layer in 1:num_layers){ #This takes about a minute to run.
  for(i in 1:8){
    land_north[[land_layer]][,i] = landscape(pd_north[[i]],K=land_layer, tseq=tseq)
    land_south[[land_layer]][,i] = landscape(pd_south[[i]],K=land_layer, tseq=tseq)
  }}


################################### Visualize landscape functions
### Visualize the landscape functions
mean_land_north = matrix(NA, nrow=length(tseq), ncol=num_layers)
sd_land_north = matrix(NA, nrow=length(tseq), ncol=num_layers)
mean_land_south = matrix(NA, nrow=length(tseq), ncol=num_layers)
sd_land_south = matrix(NA, nrow=length(tseq), ncol=num_layers)

for(land_layer in 1:num_layers){
  # Computes pointwise mean and standard deviations of north & south layers
  mean_land_north[,land_layer] = apply(land_north[[land_layer]],1, mean)
  sd_land_north[,land_layer] = apply(land_north[[land_layer]],1, sd)
  mean_land_south[,land_layer] = apply(land_south[[land_layer]],1, mean)
  sd_land_south[,land_layer] = apply(land_south[[land_layer]],1, sd)
  
  # Plots the north & south layers
  land_tib = tibble(tseq = tseq, 
                    mean_land_north = mean_land_north[,land_layer],
                    sd_land_north = sd_land_north[,land_layer],
                    mean_land_south = mean_land_south[,land_layer],
                    sd_land_south = sd_land_south[,land_layer])
  
  pp = ggplot(land_tib, aes(x=tseq)) +
    geom_ribbon(aes(ymin = mean_land_north-sd_land_north,
                    ymax = mean_land_north+sd_land_north), 
                fill = "blue", alpha = 0.5) +
    geom_ribbon(aes(ymin = mean_land_south-sd_land_south,
                    ymax = mean_land_south+sd_land_south), 
                fill = "orange", alpha = 0.5) +
    geom_line(aes(y = mean_land_north, color="North"), linewidth=2) +
    geom_line(aes(y = mean_land_south, color="South"), linewidth=2) +
    labs(x="T",y="Mean landscape", title=str_c("Layer: ", land_layer)) +
    theme_minimal() +
    scale_color_manual(values = c("blue","orange"),labels = c("North","South")) +
    theme(text = element_text(size=30),
          legend.title = element_blank(),
          legend.position = "top") 
  print(pp)
}

################################### Bootstrapped landscape function confidence bands
# Variable width bootstrap 
set.seed(98765)
B = 50 # We use B = 5000 in paper
n_samples = 8
B_quantile = 0.9 # Quantile for bootstrap confidence band
landN_boot_var = c()
landS_boot_var = c()

# Concatenate layers of landscape functions (to make one long vector per sample)
land_north_expanded = do.call(rbind, land_north) 
mean_land_north_expanded = apply(land_north_expanded,1,mean) #North pointwise mean landscape
sd_land_north_expanded = apply(land_north_expanded,1,sd) #North pointwise sd landscape

land_south_expanded = do.call(rbind, land_south) 
mean_land_south_expanded = apply(land_south_expanded,1,mean) #South pointwise mean landscape
sd_land_south_expanded = apply(land_south_expanded,1,sd) #Sout pointwise sd landscape

# Expand landscape domain for concatenated layers
tseq_expanded = c(tseq, tseq+max(tseq), tseq+2*max(tseq), tseq+3*max(tseq))

for(b in 1:B){
  sampN = sample(1:n_samples, n_samples, replace=TRUE)
  new_north = land_north_expanded[,sampN]
  new_north_mean = apply(new_north,1,mean)
  temp = which(mean_land_north_expanded == new_north_mean)
  if(length(temp)==length(tseq_expanded)){
    print("north")
    landN_boot_var[b] = 0
  }else{landN_boot_var[b] = max(((abs(mean_land_north_expanded-new_north_mean)/sd_land_north_expanded)[-temp]))}
  
  sampS = sample(1:n_samples, n_samples, replace=TRUE)
  new_south = land_south_expanded[,sampS]
  new_south_mean = apply(new_south,1,mean)
  temp2 = which(mean_land_south_expanded == new_south_mean)
  if(length(temp2)==length(tseq_expanded)){
    print("south")
    landS_boot_var[b] = 0
  }else{landS_boot_var[b] = max(((abs(mean_land_south_expanded-new_south_mean)/sd_land_south_expanded))[-temp2])}
}

# Compute quantiles of bootstrap samples for each landscape function layer
qNvar = quantile(landN_boot_var,B_quantile)
qSvar = quantile(landS_boot_var,B_quantile)


################################### Visualize variable-width bootstrap confidence bands

# Set data frame 
df1b = tibble(Tseq=tseq_expanded, mean_north = mean_land_north_expanded, mean_south = mean_land_south_expanded,
              bandN = qNvar*sd_land_north_expanded,
              bandS = qSvar*sd_land_south_expanded)

# Note:  the domain restarts for each layer
ggplot(df1b, aes(x=Tseq)) +
  geom_ribbon(aes(ymin = mean_north-bandN, ymax = mean_north+bandN), alpha = 0.2, fill = "orange") +
  geom_ribbon(aes(ymin = mean_south-bandS, ymax = mean_south+bandS), alpha = 0.2, fill = "blue") +
  geom_line(aes(y=mean_north, 
                color="North",
                linetype="North"),
            linewidth=2) +
  geom_line(aes(y=mean_south, 
                color="South", 
                linetype="South"),
            linewidth=2) +
  geom_vline(xintercept=tseq_expanded[c(1,1001,2001,3001)],
             color="brown",linetype="dashed",linewidth=1)+
  annotate("text", x = tseq_expanded[c(500,1500,2500,3500)], y = rep(.44,4), 
           label = c("k=1","k=2","k=3","k=4"),size=12,
           color="brown") +
  scale_color_manual(values=c("orange","blue"),
                     labels=c("North","South")) +
  scale_linetype_manual(values=c(1,2),
                        labels=c("North","South")) +
  guides(color = guide_legend("Group"),
         linetype = guide_legend("Group")) +
  ylab("Average landscape function") + xlab(expression((Birth+Death)/2))+
  theme_bw() +
  theme(text = element_text(size=40),
        #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top",
        legend.key.width = unit(4, "cm"),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = tseq_expanded[c(1,500, 1000,1500,2000,2500,3000,3500,4000)],
                     labels = round(c(min(tseq), rep(c(median(tseq), max(tseq)),4)),2))


###################################
