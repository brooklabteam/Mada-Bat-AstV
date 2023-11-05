library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(devtools)
library(dplyr)
library(stringr)

labs <- data.frame(
  lat = c(-19.729925, -21.118694, -17.929380),
  long = c(46.524117, 55.537113, 35.980488),
  names = c("Madagascar", "Reunion Island", "Mozambique"),
  stringsAsFactors = FALSE
)  


gg1 <- ggplot() +
  geom_polygon(data = map_data('world'),
               aes(x=long, y=lat, group=group),
               color='black', fill = "darkolivegreen4") +
  theme_bw()+
  coord_map() +
  coord_fixed(1.3,
              xlim=c(30, 56),
              ylim=c(-27, -10)) 
gg1  





gg2 <- ggplot() +
  geom_polygon(data = map_data('world'), 
       aes(x=long, y=lat, group=group),
       color='black', fill = "white") + 
  coord_map() +
  coord_fixed(1.3,
              xlim=c())
gg2
