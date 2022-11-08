library(raster)
library(rgdal)
library(tidyverse)

MAT <- raster::raster("chelsa_1970_2013/bio01.tif")
MAP <- raster::raster("chelsa_1970_2013/bio12.tif")

species_new <- read_csv("Species_parameters_template.csv")
lat_lon <- species_new %>% dplyr::select(Longitude,Latitude) %>% na.omit()
centroid_spdf <-SpatialPointsDataFrame(lat_lon,lat_lon, proj4string=MAT@crs)
MAT_values <- raster::extract(MAT,centroid_spdf,fun=mean, method="bilinear")
MAP_values <- raster::extract(MAP,centroid_spdf,fun=mean, method="bilinear")
