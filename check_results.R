#####################################################################################
############################### RESULTS #############################################
#####################################################################################

results.folder <- "./"

model1.file <- "inla_model1_results_MoliniaCaerulea_2020-07-29.RData"

# Load data
load(paste0(results.folder,model1.file))
model1 <- res.mark

### Prepare projectors ###
# Extract indices for spatial predictions 
index.pred <-
  inla.stack.index(join.stack,"pred")$data

# Mean and mode of linear predictor at mesh nodes
linpred.mean <- res.mark$summary.linear.predictor[index.pred,"mean"]

# Oxfordshire as polygon
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(coords)),'0')))
library(rgeos)
# This gets rid of TopologyException: Input geom 1 is invalid: 
# Self-intersection at or near point 502711.9115817704 174950.92257186017 
# at 502711.9115817704 174950.92257186017
# see https://gis.stackexchange.com/questions/163445/
# getting-topologyexception-input-geom-1-is-invalid-which-is-due-to-self-intersec
domainSP <- gBuffer(domainSP, byid=TRUE, width=0)

# r0 is auxiliary variable for setting the "aspect ratio" for the lattice dimensions
r0 <- diff(range(coords[, 1])) / diff(range(coords[, 2]))
# Create a mesh projector onto a grid
prj <- inla.mesh.projector(mesh.s, 
                           xlim = range(coords[, 1]),
                           ylim = range(coords[, 2]), 
                           dims = c(100, 100 / r0))

# Which points of the projector belong to the Oxfordshire domain
ov <- over(SpatialPoints(prj$lattice$loc), domainSP)

# Project mean of linear predictor field from the mesh onto a grid
prj.linpred.mean <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,linpred.mean[1:spde$n.spde + (j - 1) * spde$n.spde])
  r[is.na(ov)] <- NA
  return(r) 
})
#####

# Function to compute mark 1 prob point estimate
prob1.point.est <- function(linpred) {
  res <- exp(linpred) / (1+exp(linpred))
  return(res)
}

time.comparison = 7 # Set 7th (=2018) time knot as reference level
# Compute probability fields and ratios of probabilities
compute_prob_ratios <- function() {
  prob.mark1.linpred.mean <- vector(mode = "list", length = mesh.t$n)
  prob.mark0.linpred.mean <- vector(mode = "list", length = mesh.t$n)
  prob.ratio.linpred.mean <- vector(mode = "list", length = mesh.t$n)
  prob_ratio_means.linpred.mean <- rep(NA,mesh.t$n)
  
  for (time in 1:mesh.t$n) {
    prob.mark1.linpred.mean[[time]] <- prob1.point.est(prj.linpred.mean[[time]])
    prob.mark0.linpred.mean[[time]] <- 1 - prob.mark1.linpred.mean[[time]]
    prob.1.ratio.linpred.mean <- prob.mark1.linpred.mean[[time]] / 
      prob1.point.est(prj.linpred.mean[[time.comparison]])
    prob.0.ratio.linpred.mean <- prob.mark0.linpred.mean[[time]] / 
      (1 - prob1.point.est(prj.linpred.mean[[time.comparison]]))
    prob.ratio.linpred.mean[[time]] <- prob.1.ratio.linpred.mean / 
      prob.0.ratio.linpred.mean
    prob_ratio_means.linpred.mean[time] <- 
      mean(prob.ratio.linpred.mean[[time]],na.rm = TRUE)
  }
  result <- list()
  result$prob.ratio.linpred.mean <- prob.ratio.linpred.mean
  result$prob_ratio_means.linpred.mean <- prob_ratio_means.linpred.mean
  result$prob.mark1.linpred.mean <- prob.mark1.linpred.mean
  result$prob.mark0.linpred.mean <- prob.mark0.linpred.mean
  
  return(result)
}
prop.ratios.model1 <- compute_prob_ratios()

# Which rare and background observations belong to which time period
igr.rare.all_covs <- apply(abs(outer(data.rare.all_covs$RecYear, 
                                     mesh.t$loc, '-')), 1, which.min)
igr.bg.all_covs <- apply(abs(outer(bg_data.habitats$RecYear, 
                                   mesh.t$loc, '-')), 1, which.min)

# Plot \Xi(s,t) = p1(s,t)/p1(s,t') / p0(s,t)/p0(s,t')
for (time in 1:mesh.t$n) {
  prob.ratio <- prop.ratios.model1$prob.ratio.linpred.mean[[time]]
  main.title = paste0(time_knots[time]-tdiff+1,"-",time_knots[time]+tdiff)
  if (time == 1)
    main.title = paste0("until ",time_knots[time]+tdiff)
  if (time == 7)
    main.title = "2015-2020"
  image(x = prj$x, y = prj$y, z = prob.ratio, asp = 1, 
        xlab = '',ylab = '', zlim = range(prob.ratio,na.rm=TRUE), 
        axes = FALSE, col = cust.colors,
        main = paste0("Covered time period: ", main.title),cex=2.5,cex.main=2.5)
  lines(coords, col = "black", lwd = 1)
  image.plot(legend.only = TRUE, zlim = range(prob.ratio,na.rm=TRUE), 
             col = cust.colors,axis.args=list(cex.axis=2.5))
  points(data.rare.all_covs[igr.rare.all_covs == time, 
                            c("EastingKm","NorthingKm")], 
         pch = 19,cex=0.8)
  points(data.rare.all_covs[igr.rare.all_covs == time.comparison, 
                            c("EastingKm","NorthingKm")], 
         pch = 1,cex=0.8,col="red")
}

# Plot spatial mean \bar{Xi} = mean_s (\Xi(s,t))
last_time <- 9
ylim <- c(0.95,1.06)
plot(time_knots[1:last_time],
     na.omit(prop.ratios.model1$prob_ratio_means.linpred.mean[1:last_time]),
     pch=19,
     xlab="Time period",
     ylab=expression(bar(Xi)),
     xaxt = "n",
     ylim=ylim,
     cex=1.2,
     cex.axis=1.62,
     cex.lab = 1.65,col="blue")
abline(h=1.0,lty="dashed")
v1 <- time_knots[1:last_time]
v2 <- c("Until 1974",
        "'75-'82",
        "'83-'90",
        "'91-'98",
        "'99-'06",
        "'07-'14",
        "'15-'20",
        "2026",
        "2034")
axis(side = 1, 
     at = v1, 
     labels = v2,cex.axis=1.6)
lines(time_knots[1:7],na.omit(prop.ratios.model1$prob_ratio_means.linpred.mean[1:7]),
      lwd=2,col="blue")
lines(time_knots[7:last_time],
      na.omit(prop.ratios.model1$prob_ratio_means.linpred.mean[7:last_time]),
      lwd=2,col="blue",lty="longdash")
