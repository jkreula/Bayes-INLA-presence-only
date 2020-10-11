##################################################################################
################################# EDA ############################################
################################################################################## 
# Load libraries
library(dplyr)
library(RColorBrewer)
library(fields) # for image.plot
library(ggplot2)
# Load species data with covariates
load("data_plant_locs_all_covariates.Rda")

# How many species?
length(unique(data_plant_locs.habitats$TaxonName))

# Which years covered?
years <- sort(unique(data_plant_locs.habitats$RecYear))

# Counts per year
counts <- table(data_plant_locs.habitats$RecYear)
counts <- as.vector(counts)

# Plot plants records per year
folder <- "./Figures/"
plot.title <- paste0(folder,"plant_records_per_year_3.pdf")
pdf(file = plot.title,
    width = 7, 
    height = 4.5, 
    family = "Palatino")
par(mar = c(5,5,3,1))
plot(years,counts,
     ylim=c(0,35000),
     type='l',
     xlab="Year",
     ylab="Record count",cex=2,cex.axis=2,cex.lab=2)
dev.off()

# Plot species records per year
df.year_species <- data_plant_locs.habitats[,c("RecYear","TaxonName")]
df.year_species_agg <- df.year_species %>% 
                       group_by(RecYear) %>% 
                       summarize(n_species= n_distinct(TaxonName))
folder <- "./Figures/"
plot.title <- paste0(folder,"species_recorded_per_year.pdf")
pdf(file = plot.title,
    width = 8, 
    height = 5.5, 
    family = "Palatino")
par(mar = c(5,5,3,1))
with(df.year_species_agg,
     barplot(n_species ~ RecYear,
             ylim=c(0,1500),
             xlab="Year",
             ylab="Number of species recorded",
             cex=1.5,cex.axis=1.5,cex.lab=1.5))
dev.off()
sum(df.year_species$RecYear == 1794)

####################### RARE PLANTS ##################################
# Mask for rare plants
rare.plant.mask <- unlist(lapply(data_plant_locs.habitats$StatusOth,
                                 function(x) grepl("Oxon", x, fixed = TRUE)))
# Extract rare plants
rare.plants.df <- data_plant_locs.habitats[rare.plant.mask,]

rare.counts_per_year <- table(rare.plants.df$RecYear)
rare.counts_per_year
barplot(rare.counts_per_year)

folder <- "./Figures/"
plot.title <- paste0(folder,"rare_plant_records_per_year.pdf")
pdf(file = plot.title,
    width = 8, 
    height = 5.5, 
    family = "Palatino")
par(mar = c(5,5,3,1))
barplot(rare.counts_per_year,
        ylim=c(0,350),
        xlab="Year",
        ylab="Key species record count",cex=1.5,cex.axis=1.5,cex.lab=1.5)
dev.off()


df.rare.year_species <- rare.plants.df[,c("RecYear","TaxonName")]
df.rare.year_species_agg <- df.rare.year_species %>% 
                            group_by(RecYear) %>%
                            summarize(n_species= n_distinct(TaxonName))
folder <- "./Figures/"
plot.title <- paste0(folder,"rare_species_recorded_per_year_4.pdf")
pdf(file = plot.title,
    width = 8, 
    height = 5.5, 
    family = "Palatino")
par(mar = c(5,5,3,1))
with(df.rare.year_species_agg,
     barplot(n_species ~ RecYear,
             ylim=c(0,120),
             #xlim=c(1794,2019),
             xlab="Year",
             ylab="Number of key species recorded",
             cex=1.5,cex.axis=1.5,cex.lab=1.5))
dev.off()

### Molinia caerulea
data.rare.habitats.mc <- 
    data_plant_locs.habitats[data_plant_locs.habitats$TaxonName 
                             == "Molinia caerulea",]
data.rare.habitats.mc$CommonName[1]
folder <- "./Figures/"
file.name <- paste0(folder,"rare_species_locs_habitats_mc_2.pdf")
if (file.exists(file.name)) stop(file.name, " already exists")
pdf(file = file.name,
    width = 9, 
    height = 7, 
    family = "Palatino")
ggplot(Oxon) + geom_sf() + theme_bw() + coord_sf(datum = NA) +
    geom_point(data = data.rare.habitats.mc, 
               aes(x = Easting, y = Northing),
               size=0.8,colour="red") +
    facet_wrap(~RecYear, ncol = 10) +
    labs(x = "", y = "") +
    theme(text = element_text(size=18),
          axis.text = element_text(color="black",size=18))
dev.off()

### Iberis amara
data.rare.habitats.ia <- 
    data_plant_locs.habitats[data_plant_locs.habitats$TaxonName 
                             == "Iberis amara",]
data.rare.habitats.ia$CommonName[1]
folder <- "./Figures/"
file.name <- paste0(folder,"rare_species_locs_habitats_ia_2.pdf")
if (file.exists(file.name)) stop(file.name, " already exists")
pdf(file = file.name,
    width = 9, 
    height = 6, 
    family = "Palatino")
par(mar = c(0.1, 2.1, 4.1, 2.1))
ggplot(Oxon) + geom_sf() + theme_bw() + coord_sf(datum = NA) +
    geom_point(data = data.rare.habitats.ia, aes(x = Easting, y = Northing),
               size=0.8,colour="red") +
    facet_wrap(~RecYear, ncol = 10) +
    labs(x = "", y = "") +
    theme(text = element_text(size=18),
          axis.text = element_text(color="black",size=18))
dev.off()

##################### BACKGROUND (BASELINE) SPECIES ######################
# Choose background species
bg1 <- "Creeping Buttercup"
bg2 <- "Yorkshire-fog"
bg3 <- "Cock's-foot"
bg4 <- "False Oat-grass"
bg5 <- "Cock's-foot"
bg6 <- "Rough Meadow-grass"
bg7 <- "Hawthorn"

bg_species <- c(bg1,
                bg2,
                bg3,
                bg4,
                bg5,
                bg6,
                bg7)

# Extract background species
bg_data.habitats <- 
    data_plant_locs.habitats[data_plant_locs.habitats$CommonName %in% bg_species,]

# Plot records per year
pdf(file = "Figures/Seven_baseline_species_locs_per_year_2.pdf",
    width = 7.5, 
    height = 5, 
    family = "Palatino")
par(mar = c(0.1, 2.1, 4.1, 2.1))
ggplot(Oxon) + geom_sf() + theme_bw() + coord_sf(datum = NA) +
    geom_point(data = bg_data.habitats, aes(x = Easting, y = Northing),
               size=0.2,colour="red") +
    facet_wrap(~RecYear,ncol=14) +
    labs(x = "", y = "") +
    theme(text = element_text(size=16),
          axis.text = element_text(color="black",size=16))
dev.off()

(bg.counts.year <- as.vector(table(bg_data.habitats$RecYear)))
pdf(file = paste0("Figures/Seven_baseline_species_records_per_year.pdf"),
    width = 11,
    height = 6,
    family = "Palatino")
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1))
plot(bg.counts.year,
     ylim=c(0,1700),
     xlab="Year",
     ylab="Record count",
     cex=1.2,
     cex.axis=1.62,
     cex.lab = 1.65,
     xaxt="n",
     lwd=2)
lines(bg.counts.year,lwd=2)
axis(side = 1, 
     at = 1:length(bg.counts.year),
     labels = sort(unique(bg_data.habitats$RecYear)),cex.axis=1.6)
dev.off()

# Grouping in time
time_knots <- sort(seq(from = 1970, to = max(df.plant_data$RecYear), by = 8))
time_groups - apply(abs(outer(bg_data.habitats$RecYear, time_knots, '-')), 
                    1, which.min)
bg_data.habitats.grouped <- bg_data.habitats
bg_data.habitats.grouped$index <- time_groups
(bg.counts.grouped <- as.vector(table(bg_data.habitats.grouped$tindex)))

pdf(file = paste0("Figures/Seven_baseline_species_records_per_time_period.pdf"),
    width = 11,
    height = 6,
    family = "Palatino")
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1))
plot(bg.counts.grouped,
     ylim=c(0,6500),
     xlab="Time period",
     ylab="Record count",
     cex=1.2,
     cex.axis=1.62,
     cex.lab = 1.65,
     xaxt="n",
     lwd=2)
lines(bg.counts.grouped,lwd=2)
axis(side = 1, 
     at = 1:7, 
     labels = labs,cex.axis=1.6)
dev.off()