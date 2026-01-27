#########################################################
# Get environmental data for each sample
# to test for the effect of macroclimatic variables
# in community composition
#########################################################



######################################################
# Get a distance matrix of environmental conditions  #
######################################################
library(sf)
library(tmap)
library(raster)
library(stringr)
library(ggplot2)
library(dplyr)
library(corrplot)
library(Hmisc)


#--------------------------------------------------------------
# read sample data
#--------------------------------------------------------------
samples <- read.csv("data/sample_info_all.csv")

samples_st <- st_as_sf(samples, coords = c('longitude', 'latitude'), crs= 4326)

st_crs(samples_st) = 4326


#--------------------------------------------------------------
# Read environmental data
#--------------------------------------------------------------

# Read West Coast boundary
west.coast <- st_read(dsn= 'data/westcoast',  layer = 'west_coast')
west.coast <- st_as_sf(west.coast) #transform to sf
st_crs(west.coast) = 4326 #set projection

#californica boundary
california <- filter(west.coast, NAME == "California")

# Read bioclimatic layers
files <-  dir (path = 'data/Layers_WA')
rasters <- list()

for (i in 1:length(files)) {
  rasters[i] <- raster(paste0("data/Layers_WA/", files[i])) #loop for reading them
}

vars <- str_remove( files, '.asc')
names(rasters) = vars #naming the layers in the list

#create raster stack
environment <- stack(rasters)


#--------------------------------------------------------------
# Extract env data for collection points
#--------------------------------------------------------------

#extract environmental variables for these points
env_points <- raster::extract(environment, st_coordinates(samples_st))

#transform to df cause the output is a matrix
df_env_points <- as.data.frame(env_points)
rownames(df_env_points) <- samples$Sample


