
#####################
# Map of the samples
#####################

# useful site to check the basemaps: https://cran.r-project.org/web/packages/maptiles/maptiles.pdf

Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj")  
Sys.setenv(GDAL_DATA = "/opt/homebrew/share/gdal")

library(tmap)
library(sf)
library(maptiles)
library(viridis)
library(sf)

# some checks
sf::sf_proj_search_paths()
file.exists(file.path(Sys.getenv("PROJ_LIB"), "proj.db"))


# Read samples
samples <- read.csv("data/sample_info_all.csv", header = T)

# Make an sf object
samples_sf <- st_as_sf(samples, coords = c('longitude', 'latitude'), crs=4326)
st_crs(samples_sf)

# Make a dynamic map
tmap_mode("view")

tm_shape(samples_sf) + tm_dots(fill = "Geographic_region", size = 0.65)


# make a static map
tmap_mode("plot")

tm_basemap("Esri.WorldImagery") +
    tm_shape(samples_sf) +
    tm_dots(fill = "yellow", size = 0.7) +
    tm_compass(type="8star", position=c("left", "bottom"), show.labels = 3, text.color = "white")


pal <- c("Centralcal"="#440154FF", 
         "Desert"="#EB5760FF", 
         "Northcal"="#2E6E8EFF", 
         "Peninsular ranges"="#1F9E89FF", 
         "SCI"="#982D80FF", 
         "Sierras"="#71CF57FF", 
         "Transverse ranges"="#FDE725FF")

# 
# tm_basemap("Esri.WorldImagery") +
#     tm_shape(samples_sf) +
#     tm_dots(fill = "Geographic_region", fill_alpha= 0.9,
#         size = 0.8,
#         fill.scale = tm_scale(values = pal)) +
#     tm_compass(type = "8star",
#         position = c("left", "bottom"),
#         show.labels = 3,
#         text.color = "white")


# Map for figure 1
map <- tm_basemap("Esri.WorldTerrain") +
    tm_shape(samples_sf) +
    tm_dots(fill = "Geographic_region", fill_alpha= 0.9,
            size = 0.8,
            fill.scale = tm_scale(values = pal)) +
    tm_compass(type = "8star",
               position = c("left", "bottom"),
               show.labels = 3,
               text.color = "black") +
    tm_layout(legend.position = c("right", "top")) 

map

# save
png("plots/map.png", res = 300, units = "cm", width = 10, height = 15)
map
dev.off()




