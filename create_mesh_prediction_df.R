### © Juha Kreula

require(sf)
require(tidyr) # for drop_na
require(raster) # for handling raster objects
require(spatstat) # for nearest.pixel
require(maptools) # for as.im.RasterLayer()

# Load covariate data frame?
load.covar.df <- ifelse(file.exists("all_covariate_data.RData"),TRUE,FALSE) 

if (load.covar.df) {
  load("all_covariate_data.RData")
} else {
  # Habitats (change path if needed)
  data_habitats <- st_read("./Data/Oxon-Plant-Records/Oxon-Habitats.shp")
  
  # Topographical data (change path if needed)
  folder <- "./Data/topographical_data/"
  
  aspect.file <- paste0(folder,"aspect.tif")
  elevation.file <- paste0(folder,"elevation.tif")
  geology.file <- paste0(folder,"geology.tif")
  slope.file <- paste0(folder,"slope.tif")
  oxon.extent <- extent(c(Oxon_limits$xmin,
                          Oxon_limits$xmax,
                          Oxon_limits$ymin,
                          Oxon_limits$ymax))
  data.aspect <- raster(aspect.file)
  data.aspect.oxon <- crop(data.aspect, oxon.extent)
  data.elevation <- raster(elevation.file)
  data.elevation.oxon <- crop(data.elevation, oxon.extent)
  data.geology <- raster(geology.file)
  data.geology.oxon <- crop(data.geology, oxon.extent)
  data.slope <- raster(slope.file)
  data.slope.oxon <- crop(data.slope, oxon.extent)
  # Save for future
  save(data_habitats,
       data.aspect.oxon,
       data.elevation.oxon,
       data.geology.oxon,
       data.slope.oxon,
       file = "all_covariate_data.RData")
}


# Use mesh nodes as prediction points
pred.xy <- st_as_sf(data.frame(x = (mesh.s$loc[,1]*1000), y = (mesh.s$loc[,2])*1000),  
                    coords = c("x", "y")  )  
st_crs(pred.xy) <- st_crs(27700)
# Find which polygon (corresponding to a habitat area on map) each observation point belongs to
# and create a data frame
predxy.habit.intersection <- st_intersects(data_habitats$geometry,pred.xy)
df.predxy.habit.intersection <- as.data.frame(predxy.habit.intersection)
# Habitat index is the row number in data_habitats, 
# PointIndex is the row number in data_plant_locs
colnames(df.predxy.habit.intersection) <- c("HabitatIndex","PointIndex") 

# Placeholders for new columns or features for habitat in data_plant_locs

pred.xy$HabitatTypeLong <- rep(NA, dim(pred.xy)[1])

pred.xy$HabitatTypeLong[df.predxy.habit.intersection$PointIndex] = 
  data_habitats$IHSHABITAT[df.predxy.habit.intersection$HabitatIndex]

pred.xy$HabitatTypeShort <- sapply(pred.xy$HabitatTypeLong, 
                                   function(x) substr(x, start = 1, stop = 2))


pred.xy.HabitatInfo <- pred.xy[!is.na(pred.xy$HabitatTypeShort),]
pred.xy.NoHabitatInfo <- pred.xy[is.na(pred.xy$HabitatTypeShort),]


create_pred_df <- function(pred) {
  result <- data.frame(CommonName = NA, 
                       TaxonName = NA, 
                       RecYear = NA, 
                       Easting = st_coordinates(pred$geometry)[,1],
                       Northing = st_coordinates(pred$geometry)[,2],
                       HabitatTypeLong = pred$HabitatTypeLong,
                       HabitatTypeShort = pred$HabitatTypeShort,
                       Aspect = NA,
                       Elevation = NA,
                       Geology = NA,
                       Slope = NA,
                       IS_DATA = FALSE)
  return (result)
}


df.pred.xy.HabitatInfo <- create_pred_df(pred.xy.HabitatInfo)
df.pred.xy.NoHabitatInfo <- create_pred_df(pred.xy.NoHabitatInfo)

df.plant_data$IS_DATA <- TRUE
df.all.HabitatInfo <- rbind(df.plant_data,df.pred.xy.HabitatInfo)

euclidean_dist <- function(x1,y1,x2,y2) {
  dist <- sqrt( (x1 - x2)^2 + (y1 - y2)^2  )
  return(dist)
}

find.imputed.habitats <- function(nohabitat.locs,df.habitats)
{
  imputed.habitats <- apply(nohabitat.locs, 1, function(row) {
    row.dist <- euclidean_dist(row[1],
                               row[2],
                               df.habitats$Easting,
                               df.habitats$Northing)
    arg.min.dist <- which.min(row.dist)
    habitat <- df.habitats$HabitatTypeShort[arg.min.dist]
    return(habitat)
  })
  return(imputed.habitats)
}

nohabitat.locs <- cbind(df.pred.xy.NoHabitatInfo$Easting,
                        df.pred.xy.NoHabitatInfo$Northing)


imputed.habitats <- find.imputed.habitats(nohabitat.locs,
                                          df.all.HabitatInfo)

df.pred.xy.ImputedHabitatInfo <- df.pred.xy.NoHabitatInfo
df.pred.xy.ImputedHabitatInfo$HabitatTypeShort <- imputed.habitats

df.pred.xy <- rbind(df.pred.xy.HabitatInfo,df.pred.xy.ImputedHabitatInfo)

head(df.pred.xy)
n.preds <- dim(df.pred.xy)[1]

# Include continuous covariates to data frame
extract.covariate.vals <- function(df,X) {
  data.im <- maptools::as.im.RasterLayer(X)
  x <- df$Easting
  y <- df$Northing
  pixels <- nearest.pixel(x,y,data.im)
  vals <- sapply(1:length(x),function(i) {X[pixels$row[i],pixels$col[i]] } )
  return(vals)
}

df.pred.xy$Aspect <- extract.covariate.vals(df.pred.xy,data.aspect.oxon)
df.pred.xy$Elevation <- extract.covariate.vals(df.pred.xy,data.elevation.oxon)
df.pred.xy$Geology <- extract.covariate.vals(df.pred.xy,data.geology.oxon)
df.pred.xy$Slope <- extract.covariate.vals(df.pred.xy,data.slope.oxon)

df.pred.xy$Geology <- as.factor(df.pred.xy$Geology)