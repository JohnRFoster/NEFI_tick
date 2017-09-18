library(rgbif)
library(rinat)
library(tidyverse)
library(maps)

# White-footed mouse = P. leucopus
# Eastern chipmunk = T. striatus
# Short-tailed shrew = B. brevicauda
# Masked shrew = S. cinereus
# Black-legged tick = I. scapularis
# Western black-legged tick = I. pacificus

spp_list <- c("Peromyscus leucopus", "Tamias striatus", "Blarina brevicauda", "Sorex cinereus",
              "Ixodes scapularis", "Ixodes pacificus")

#---------------------------------------------
# GBIF Data Intake

# data fields to be returned in occ_search not specified in occ_search call
fields <- c("scientificName", "eventDate", "year", "month", "decimalLatitude", "decimalLongitude",
            "elevation", "protocol", "publishingOrg", "issue")

# GBIF Search function
gbif_occ_func <- function(species){
  occ_search(scientificName = species, hasCoordinate = TRUE, country = "US", fields = fields,
             basisOfRecord = "HUMAN_OBSERVATION", return = "data", limit = 5000)
}

gbif_occ <- sapply(spp_list, gbif_occ_func) %>% t()

#-------------------------------------------------
# iNat Data Intake

# data fileds to be returned in get_inat_obs
inat_return <- c("scientific_name", "datetime", "latitude", 
                 "longitude", "id", "observed_on", "positional_accuracy")

# iNat Search function
inat_occ_func <- function(species){
  get_inat_obs(taxon_name = species, geo = TRUE, quality = "research")
}

inat_occ <- sapply(spp_list, inat_occ_func) %>% t() 
inat_occ <- inat_occ[, inat_return]

#--------------------------------------------------