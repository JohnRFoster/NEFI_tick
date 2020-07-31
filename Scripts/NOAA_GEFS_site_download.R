# need to set this for CRON, will default to geo home instead of projectnb
setwd("/projectnb/dietzelab/fosterj/NEFI_tick")

library(rnoaa)
library(ncdf4)
library(lubridate)
library(plantecophys)
library(tidyr)
library(dplyr)
library(PEcAn.logger)
library(PEcAn.remote)

source("Functions/NOAA_GEFS.R")
overwrite <- FALSE

cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")

today <- Sys.Date()

# define sitenames
sitename <- c(
  "HarvardForest",
  "CaryInstitute"
)

# site level storage
outfolder <- c(
  "../GEFS/HarvardForest",
  "../GEFS/Cary"
)

# latitude for each site
lat.in <- c(
  42.5369,     # Harvard Forest
  41.7851      # Cary Institute 
)

# longitude for each site
lon.in <- c(
  -72.17266,   # Harvard Forest
  -73.7338     # Cary Institute
)

# timezone for each site
tz <- c(
  "America/New_York",  # Harvard Forest
  "America/New_York"   # Cary Institute
)

# download GEFS for each site and save results
for(site in seq_along(sitename)){
  
  # don't need to run if already downloaded
  check.folder <- paste("NOAA_GEFS", sitename[site], today, sep = ".")
  check <- file.path(outfolder[site], check.folder)
  if(dir.exists(check) & !overwrite){
    cat(sitename[site], "already downloaded!\n")
    next
  }
  
  # get GEFS
  cat("===== Attempting", sitename[site], "download =====\n")
  results <- NOAA_GEFS(
    outfolder[site],
    lat.in[site],
    lon.in[site],
    sitename[site],
    time.zone = tz[site],
  )
  
  save(results, 
       file = paste0(outfolder[site], "/results.RData"))
  
  cat("=====", sitename[site], "downloaded =====\n")
  
}

cat("\n ----- END ----- \n\n")

