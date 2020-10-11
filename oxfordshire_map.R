#Clear workspace
rm(list=ls())
# Load libraries
library(lwgeom)  
library(raster) 
library(sf)  
library(dplyr) 
library(ggplot2)

# Fetch map for Great Britain
m <- raster::getData(name = "GADM", country = "GB", level = 2)

# Extract Oxfordshire
Oxon_map = m[m$NAME_2 == "Oxfordshire", ]

# Transform object to POLYGON
Oxon_pol <- Oxon_map %>%
  st_as_sf() %>%
  st_cast("POLYGON")

# Plot map
pdf(file = "oxfordshire_map_longlat.pdf",
    width = 5,
    height = 5,
    family = "Palatino")
par(mfrow = c(1, 1), mar = c(1, 1, 2, 3))
ggplot(Oxon_pol) + geom_sf() + theme_bw() +
  xlab("Longitude")+
  ylab("Latitude") +
  theme(text = element_text(size=18),
        axis.text = element_text(color="black",size=18))
dev.off()

# Transform map to UTM coordinates (i.e. easting and northing instead of longitude and latitude)
# The EPSG code for Great Britain is 27700
Oxon <- Oxon_pol %>% st_transform(27700)
oxon_boundary_coords <- st_coordinates(Oxon)[,1:2]
Oxon_limits <- st_bbox(Oxon)

# Save UTM Oxon object
save(Oxon,oxon_boundary_coords,Oxon_limits,file = "Oxon_map.RData")

# Plot Oxfordshire map and show the location of the city of Oxford
pdf(file = "oxfordshire_map_UTM.pdf",
   width = 5.1,
   height = 5,
   family = "Palatino")
print(ggplot(Oxon) + geom_sf() + theme_bw() + coord_sf(datum = st_crs(Oxon)) +
  geom_point(aes(x = 451342, y = 206181),size=2) +
  geom_label(label="Oxford",aes(x = 451342, y = 201181),size=7) +
  scale_x_continuous(breaks=c(420000,440000,460000,480000)) +
  xlab("Easting (m)")+
  ylab("Northing (m)") +
  theme(text = element_text(size=22),
        axis.text = element_text(color="black",size=22)) +
  theme(plot.margin=unit(c(0.1,1,0,0),"cm")))
dev.off()