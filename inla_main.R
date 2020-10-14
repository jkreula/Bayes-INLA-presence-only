### © Juha Kreula

### Load libraries
library(raster)
library(sf)
library(lwgeom)  
library(dplyr) 
library(ggplot2)
library(RColorBrewer) # For plotting with custom colors
library(spatstat) # for nearest.pixel
library(maptools) # for as.im.RasterLayer()
library(readr) # for extract.numeric
library(fields) # for image.plot
library(tidyr) # for drop_na
library(raster) # for handling raster objects
library(spatstat) # for nearest.pixel
library(maptools) # for as.im.RasterLayer()
library(fields) # for image.plot


########## Data load and focus species ######################

# Load Oxon map
load("Oxon_map.RData")

# Load species data with covariates
load("data_plant_locs_all_covariates.Rda")

# Choose columns to keep in data
df.plant_data <- data_plant_locs.habitats[,c("CommonName",
                                             "TaxonName",
                                             "RecYear",
                                             "Easting",
                                             "Northing",
                                             "HabitatTypeLong",
                                             "HabitatTypeShort",
                                             "Aspect",
                                             "Elevation",
                                             "Geology",
                                             "Slope")]



# Choose rare species to be investigated. Use Latin name.
rare.investigate <- "Molinia caerulea"
stopifnot(rare.investigate %in% unique(df.plant_data$TaxonName))



########### Background species #######################

# Choose background species
bg1 <- "Creeping Buttercup"
bg2 <- "Yorkshire-fog"
bg3 <- "Cock's-foot"
bg4 <- "False Oat-grass"
bg5 <- "Cock's-foot"
bg6 <- "Rough Meadow-grass"
bg7 <- "Hawthorn"

# Vector of background species
bg_species <- c(bg1,
                bg2,
                bg3,
                bg4,
                bg5,
                bg6,
                bg7
)

# Data frame of background species
bg_data.habitats <- df.plant_data[df.plant_data$CommonName %in% bg_species,]

# Add mark to bg point data
bg_data.habitats$Mark <- 0

# Extract dataframe for rare species
data.rare.all_covs <- df.plant_data[df.plant_data$TaxonName == rare.investigate,]
data.rare.all_covs$Mark <- 1

# Combine background and key species data frames
data.bg_and_rare <- rbind(data.rare.all_covs,bg_data.habitats)

# Find species observation locations for all times in kilometres
data.bg_and_rare$EastingKm <- data.bg_and_rare$Easting/1000
data.bg_and_rare$NorthingKm <- data.bg_and_rare$Northing/1000
sp_locs <- cbind(data.bg_and_rare$EastingKm, data.bg_and_rare$NorthingKm)

###################################################################################
################################## INLA  ##########################################
###################################################################################
if(!require(INLA)) {
  install.packages("INLA", 
                   repos=c(getOption("repos"), 
                           INLA="https://inla.r-inla-download.org/R/testing"), 
                   dep=TRUE)
  library(INLA)
}

#Needed on grey01.cpu
inla.setOption(mkl=TRUE)

# Function to create an extended triangular mesh around Oxfordshire
create_spatial_mesh <- function(domain_boundary,locations,offset,max_edge,cutoff) {
  # Change Oxfordshire coordinates from m to km
  coords <- domain_boundary / 1000
  
  # Create a non-convex hull around Oxfordshire
  bnd <- inla.nonconvex.hull(coords)
  # Create a spatial triangular mesh for finite element method
  mesh.s <- inla.mesh.2d(loc = locations,
                         boundary = bnd,
                         offset = c(2,5),
                         # Max triangle edge length in the interior is 5km, 
                         # and in the exterior 20km
                         max.edge = c(5, 20), 
                         # smallest triangle edge length is 1km
                         cutoff = 1)
  return(mesh.s)
}

# Spatial mesh
mesh.s <- create_spatial_mesh(oxon_boundary_coords,sp_locs,offset=c(2,5),max_edge=c(5, 20),cutoff=1)

# Function to create 1D temporal mesh (use coarser time grid for computational efficiency)
create_temporal_mesh <- function(start,finish,time_step,num_forecast_points=0) {
  # Add forecast points to create last time knot
  last <- finish + num_forecast_points*time_step
  
  time_knots <- sort(seq(from = start, 
                         to = last,
                         by = time_step))
  # Create mesh
  mesh.t <- inla.mesh.1d(time_knots)
  
  return(mesh.t)
}

# Create a 1D temporal mesh
mesh.t <- create_temporal_mesh(start=1970,
                               finish=max(df.plant_data$RecYear),
                               time_step=8,
                               num_forecast_points=2)
# Number of time knots
(k <- mesh.t$n)

#################### Create stacks ##############################

# Estimation projector matrix
Ast <- inla.spde.make.A(mesh = mesh.s,
                        loc = sp_locs,
                        n.group = mesh.t$n,
                        group = data.bg_and_rare$RecYear,
                        group.mesh = mesh.t)

# Create space-time index
index <- inla.spde.make.index("s", spde$n.spde, n.group = k)

# Mark response variable
mark <- data.bg_and_rare$Mark

# Extract covariate values at observation locations
Habs.sp <- as.factor(data.bg_and_rare$HabitatTypeShort)
Aspect.sp <- as.numeric(data.bg_and_rare$Aspect)
Elevation.sp <- as.numeric(data.bg_and_rare$Elevation)
Geology.sp <- as.factor(data.bg_and_rare$Geology)
Slope.sp <- as.numeric(data.bg_and_rare$Slope)

# Function to standardise covariates
standardise_covariates <- function(x, mean = NULL, sd = NULL) {
  if (is.null(mean) || is.null(sd)) {
    x.std <- (x - mean(x)) / sd(x)
  } else {
    x.std <- (x - mean)/sd
  }
  return(x.std)
}

# Standardise continuous covariates
Aspect.sp.std <- standardise_covariates(Aspect.sp)
Elevation.sp.std <- standardise_covariates(Elevation.sp)
Slope.sp.std <- standardise_covariates(Slope.sp)

# Extract means and stds for standardising predictions
Aspect.mean <- mean(Aspect.sp)
Elevation.mean <- mean(Elevation.sp)
Slope.mean <- mean(Slope.sp)

Aspect.sd <- sd(Aspect.sp)
Elevation.sd <- sd(Elevation.sp)
Slope.sd <- sd(Slope.sp)

########## Prediction at mesh nodes

source("create_mesh_prediction_df.R")

# Repeate covariate values over time (assumed not to have changed)
Habs.pred <- as.factor(rep(df.pred.xy$HabitatTypeShort,mesh.t$n))
Aspect.pred <- as.numeric(rep(df.pred.xy$Aspect,mesh.t$n))
Elevation.pred <- as.numeric(rep(df.pred.xy$Elevation,mesh.t$n))
Geology.pred <- as.factor(rep(df.pred.xy$Geology,mesh.t$n))
Slope.pred <- as.numeric(rep(df.pred.xy$Slope,mesh.t$n))

# Standardise at mesh nodes
Aspect.pred.std <- standardise_covariates(Aspect.pred,Aspect.mean,Aspect.sd)
Elevation.pred.std <- standardise_covariates(Elevation.pred,Elevation.mean,Elevation.sd)
Slope.pred.std <- standardise_covariates(Slope.pred,Slope.mean,Slope.sd)

# Create estimation stack
stack.est <- inla.stack(
  data = list(mark = mark),
  A = list(Ast,1),
  effects = list(index, list(Intercept = 1,
                             Habitat = Habs.sp, 
                             Aspect = Aspect.sp.std, 
                             Elevation = Elevation.sp.std,
                             Geology = Geology.sp, 
                             Slope = Slope.sp.std)),
  tag = "est"
)

# Projection from mesh to itself is identity
A.pred <- Diagonal(n = k*spde$n.spde)

# Create prediction stack 
stack.pred <- inla.stack(data = list(mark=NA),
                         A = list(A.pred,1),
                         effects = list(index, 
                                        list(Intercept = 1,
                                             Habitat = Habs.pred, 
                                             Aspect = Aspect.pred.std, 
                                             Elevation = Elevation.pred.std,
                                             Geology = Geology.pred, 
                                             Slope = Slope.pred.std)),
                         tag="pred")

# Join stacks together
join.stack <- inla.stack(stack.est,stack.pred)

############ Priors #############################

# PC prior for AR(1) process, base model rho = 1
pcrho <- list(prior = 'pc.cor1', param = c(0,0.9))

# Log-gamma prior for log precision for categoricals
prec.init <- 5e-5
prec.prior <- list(prec = list(prior = "loggamma", 
                               param = c(1, prec.init), 
                               fixed = FALSE)) 

# PC prior for Matern field on spatial mesh
range <- 10
range.alpha <- 0.1
sigma <- 0.35
sigma.alpha <- 0.05
spde <- inla.spde2.pcmatern(mesh = mesh.s,
                            prior.range = c(range, range.alpha),
                            prior.sigma = c(sigma, sigma.alpha))

# Model formula
formula.mark <- mark ~ 0 + Intercept + 
  f(Habitat, model = "iid", hyper = prec.prior, constr = TRUE) + 
  Aspect + 
  Elevation + 
  f(Geology, model = "iid", hyper = prec.prior, constr = TRUE) + 
  Slope + 
  f(s, model = spde, group = s.group, control.group = 
      list(model = "ar1", hyper = list(theta = pcrho)))

######### Run INLA #################
res.mark <- inla(formula.mark, 
                 family = 'binomial', 
                 data = inla.stack.data(join.stack),
                 Ntrials = 1,
                 control.predictor = list(compute = TRUE, 
                                          A = inla.stack.A(join.stack), link = 1),
                 control.inla = list(strategy = 'simplified.laplace', 
                                     int.strategy = 'ccd'),
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))

####### Save results #######

# Auxiliary function to modify taxon name
# Molinia caerulea -> MoliniaCaerulea for saving results
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
# Get taxon name in TitleCase
rare.species.taxonName <- gsub(" ", "", CapStr(rare.investigate))
date <- Sys.Date()

# Save filename
save.filename <- paste0("inla_model1_results_",rare.species.taxonName,
                        "_",date,".RData")
# Save workspace
save.image(file = save.filename)