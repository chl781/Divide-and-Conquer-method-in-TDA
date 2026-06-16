PlotLakes = function(){
part_coord = readRDS("partition_coordinates.rds")
names(part_coord)

my_lat_breaks_n <- part_coord$north_lat       
my_lon_breaks_n <- part_coord$north_long    
my_lat_breaks_s <- part_coord$south_lat  
my_lon_breaks_s <- part_coord$south_long

# Plot of lake locations with partitions
ggplot(lakes, aes(x=long, y=lat)) +
  geom_point() +
  ####North
  geom_rect(aes(xmin=min(my_lon_breaks_n), xmax=max(my_lon_breaks_n), 
                ymin=min(my_lat_breaks_n), ymax= max(my_lat_breaks_n)), 
            fill=NA, color="orange",linewidth=1.25) +
  geom_segment(aes(x=my_lon_breaks_n[2], xend=my_lon_breaks_n[2], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[3], xend=my_lon_breaks_n[3], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[4], xend=my_lon_breaks_n[4], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[5], xend=my_lon_breaks_n[5], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[6], xend=my_lon_breaks_n[6], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[7], xend=my_lon_breaks_n[7], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_n[8], xend=my_lon_breaks_n[8], 
                   y=min(my_lat_breaks_n), yend=max(my_lat_breaks_n)), 
               color="orange",linewidth=1.25)+
  ###South
  geom_rect(aes(xmin=min(my_lon_breaks_s), xmax=max(my_lon_breaks_s), 
                ymin=min(my_lat_breaks_s), ymax= max(my_lat_breaks_s)), 
            fill=NA, color="blue",linewidth=1.25) +
  geom_segment(aes(x=my_lon_breaks_s[2], xend=my_lon_breaks_s[2], 
                   y=min(my_lat_breaks_s), yend=max(my_lat_breaks_s)), 
               color="blue",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_s[3], xend=my_lon_breaks_s[3], 
                   y=min(my_lat_breaks_s), yend=max(my_lat_breaks_s)), 
               color="blue",linewidth=1.25)+
  geom_segment(aes(x=my_lon_breaks_s[4], xend=my_lon_breaks_s[4], 
                   y=min(my_lat_breaks_s), yend=max(my_lat_breaks_s)), 
               color="blue",linewidth=1.25)+
  geom_segment(aes(x=min(my_lon_breaks_s), xend=max(my_lon_breaks_s), 
                   y=my_lat_breaks_s[2], yend=my_lat_breaks_s[2]), 
               color="blue",linewidth=1.25)+
  theme_minimal() +
  labs(x="Longitude", y="Latitude") +
  theme(text = element_text(size=50)) +
  scale_x_continuous(expand = expansion(add = c(1, 1)),
                     breaks = c(-93,-91,-89,-87),
                     labels = c(-93,-91,-89,-87)
  )
}
