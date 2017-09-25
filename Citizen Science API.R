library(rgbif)
library(rinat)
library(tidyverse)
library(maps)
library(mapr)
library(sp)
library(ggmap)
library(raster)
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

# Map - Citizen Science Tick Dataset - Combined GBIF & iNat

# all tick observations from GBIF and iNat, then select 'name', 'latitude', and 'longitude' columns
tick_occ <- bind_rows(gbif_occ, inat_occ) %>% select(name, latitude, longitude) %>% 
  filter(name == "Ixodes pacificus" | name == "Ixodes scapularis") # extract just tick occurrences

# all small mammal (smam) observation from GBIF and iNat, select 'name', 'latitude', and 'longitude' columns
smam_occ <- bind_rows(gbif_occ, inat_occ) %>% select(name, latitude, longitude) %>% 
  filter(name == "Peromyscus leucopus" | name == "Tamias striatus" | name == "Blarina brevicauda" |
  name == "Sorex cinereus") # extract just smam occurrences

tick_pts <- map_ggplot(tick_occ, map = "usa", size = 2) # tick occurrences on USA map
smam_pts <- map_ggplot(smam_occ, map = "usa", size = 2) # smam occurrences on USA map


coordinates(tick_occ) <- c("longitude", "latitude")
coordinates(smam_occ) <- c("longitude", "latitude")

#tick_box <- gridlines(tick_occ, easts = pretty(bbox(tick_occ)[1, ]), norths = pretty(bbox(tick_occ)[2, ]), ndiscr = 1000)
  