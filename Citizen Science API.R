library(rgbif)
library(rinat)
library(tidyverse)
library(sp)
library(ggmap)
library(raster)
library(rasterVis)
library(measurements)

# White-footed mouse = P. leucopus
# Eastern chipmunk = T. striatus
# Short-tailed shrew = B. brevicauda
# Masked shrew = S. cinereus
# Black-legged tick = I. scapularis
# Western black-legged tick = I. pacificus

spp_list <- c("Peromyscus leucopus", "Tamias striatus", "Blarina brevicauda", "Sorex cinereus",
              "Ixodes scapularis", "Ixodes pacificus")

tick_list <- c("Ixodes scapularis", "Ixodes pacificus")
smam_list <- c("Peromyscus leucopus", "Tamias striatus", "Blarina brevicauda", "Sorex cinereus")

#---------------------------------------------

# GBIF Data Intake

# data fields to be returned in occ_search not specified in occ_search call
fields <- c("scientificName", "eventDate", "year", "month", "decimalLatitude", "decimalLongitude",
            "elevation", "protocol", "publishingOrg", "issue", "identifier", "recordedBy")

# GBIF Search function for georeferenced occurrences made by a human
gbif_occ_func <- function(species){
  occ_search(scientificName = species, hasCoordinate = TRUE, country = "US", fields = fields,
             basisOfRecord = "HUMAN_OBSERVATION", return = "data", limit = 10000)
}

# GBIF loop, applies search function to each animal in spp_list, then adds results to gbif_occ
gbif_occ <- NULL
for(i in spp_list){
  if(i == "Peromyscus leucopus"){
    gbif_occ <- gbif_occ_func("Peromyscus leucopus")
  }
  if(i == "Tamias striatus"){
    x <- gbif_occ_func("Tamias striatus")
    gbif_occ <- rbind(gbif_occ, x)
  }
  if(i == "Blarina brevicauda"){
    x <- gbif_occ_func("Blarina brevicauda")
    gbif_occ <- rbind(gbif_occ, x)
  }
  if(i == "Sorex cinereus"){
    x <- gbif_occ_func("Sorex cinereus")
    gbif_occ <- rbind(gbif_occ, x)
  }
  if(i == "Ixodes pacificus"){
    x <- gbif_occ_func("Ixodes pacificus")
    gbif_occ <- rbind(gbif_occ, x)
  }
  if(i == "Ixodes scapularis"){
    x <- gbif_occ_func("Ixodes scapularis")
    gbif_occ <- rbind(gbif_occ, x)
  }
}

gbif_occ <- gbif_occ %>% separate(scientificName, into =  c("name", "1", "2"), sep = " ") %>% 
  unite(col = "name", sep = " ", ... = c("name", "1"), remove = TRUE)

colnames(gbif_occ)[c(4, 5)] <- c("longitude", "latitude") # change the name of columns 2 - 4 to
                                                          # play nice with mapr


#-------------------------------------------------

# iNat Data Intake

# data fileds to be subset in inat_occ
inat_return <- c("name", "datetime", "latitude", 
                 "longitude", "id", "observed_on", "positional_accuracy")

# bounding box for inat query
# lat long info from https://en.wikipedia.org/wiki/List_of_extreme_points_of_the_United_States
bounds_deg <- c("24 31 15", "-124 46 18.1", "49 23 04", "-66 56 59.2") 
bounds <- conv_unit(bounds_deg, from = "deg_min_sec", to  = "dec_deg") %>% as.numeric() 


# iNat Search function for occurrences that are georeferenced and of 'research' quality
inat_occ_func <- function(species){
  get_inat_obs(taxon_name = species, geo = TRUE, quality = "research", maxresults = 10000, bounds = bounds)
}

# iNat loop, applies search function to each animal in spp_list, then adds results to inat_occ
inat_occ <- NULL
for(i in spp_list){
  if(i == "Peromyscus leucopus"){
    inat_occ <- inat_occ_func("Peromyscus leucopus")
  }
  if(i == "Tamias striatus"){
    x <- inat_occ_func("Tamias striatus")
    inat_occ <- rbind(inat_occ, x)
  }
  if(i == "Blarina brevicauda"){
    x <- inat_occ_func("Blarina brevicauda")
    inat_occ <- rbind(inat_occ, x)
  }
  if(i == "Sorex cinereus"){
    x <- inat_occ_func("Sorex cinereus")
    inat_occ <- rbind(inat_occ, x)
  }
  if(i == "Ixodes pacificus"){
    x <- inat_occ_func("Ixodes pacificus")
    inat_occ <- rbind(inat_occ, x)
  }
  if(i == "Ixodes scapularis"){
    x <- inat_occ_func("Ixodes scapularis")
    inat_occ <- rbind(inat_occ, x)
  }
} 
inat_occ <- inat_occ[!(inat_occ$id %in% gbif_occ$identifier), ]   # removes iNat data duplicated in GBIF data

inat_occ <- inat_occ %>% separate(scientific_name, into =  c("name", "1", "2"), sep = " ") %>% 
              unite(col = "name", sep = " ", remove = TRUE, ... = c("name", "1"))

inat_occ <- inat_occ[ , inat_return] # the above two calls rename and remove some columns to help the data 
                                     # play nice with mapr

#--------------------------------------------------

# Combine GBIF and iNat datasets for ticks and smams, then seperate ticks and smams into their own data frame

# all tick observations from GBIF and iNat, then select 'name', 'latitude', and 'longitude' columns
tick_occ <- bind_rows(gbif_occ, inat_occ) %>% dplyr::select(name, latitude, longitude) %>% 
  filter(name == "Ixodes pacificus" | name == "Ixodes scapularis") # extract just tick occurrences

# all small mammal (smam) observation from GBIF and iNat, select 'name', 'latitude', and 'longitude' columns
smam_occ <- bind_rows(gbif_occ, inat_occ) %>% dplyr::select(name, latitude, longitude) %>% 
  filter(name == "Peromyscus leucopus" | name == "Tamias striatus" | name == "Blarina brevicauda" |
  name == "Sorex cinereus") # extract just smam occurrences

#---------------------------------------------------

# get a map of contiguous USA boarder
usa.map <- getData(name = "GADM", country = "usa", level = 0) # USA boarder map; class "SpatialPolygonsDataFrame"
usa.map <- crop(usa.map, extent(-130, -65, 23, 52))

#---------------------------------------------------

# function 'count.raster' creates a 'RasterLayer' object from input data and counts number of occurrences in each grid cell

count.raster <- function(data, ncol, nrow) {
  data[c("name", "longitude", "latitude")]        # subset input data
  coordinates(data) <- c("longitude", "latitude") # convert data to class "SpatialPointsDataFrame"
  r <- raster()                                   # initialize raster
  extent(r) <- extent(data)                       # set extent of raster equal to input data
  ncol(r) <- ncol
  nrow(r) <- nrow
  rasterize(data, r, field = data$name, fun = 'count') # creates "RasterLayer" and counts points inside each cell
}

tick.raster <- count.raster(data = tick_occ, ncol = 40, nrow = 25) # create RasterLayer of tick occurrences with grid cells
smam.raster <- count.raster(data = smam_occ, ncol = 40, nrow = 25) # create RasterLayer of tick occurrences with grid cells

#---------------------------------------------------

# count.plot creates a "trellis" object 
count.plot <- function(raster){
 levelplot(raster) +                             # creates plot of input raster
  layer(sp.polygons(usa.map))                    # adds usa boarder to map 
}

tick.plot <- count.plot(tick.raster)  # plots tick raster with usa.map, a trellis object
smam.plot <- count.plot(smam.raster)  # plots smam raster with usa.map, a trellis object

