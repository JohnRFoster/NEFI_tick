library(ncdf4)
library(rnoaa)
library(lubridate)
library(plantecophys)
library(tidyr)
library(dplyr)
library(PEcAn.logger)
library(PEcAn.remote)

setwd("/projectnb/dietzelab/fosterj/NEFI_tick")
source("Functions/NOAA_GEFS.R")


start_date <- Sys.time() # start date/time for all sites 
cat("Start date:", format(start_date, "%Y-%m-%d %H:%M"), "\n")

#------------------------
#    Cary Institute
#------------------------
lat.in <- 41.7851
lon.in <- -73.7338
sitename <- "CaryInstitute"
outfolder <- "../GEFS/Cary"
tz <- "America/New_York"
end_date <- as.POSIXct(start_date, tz = tz) + lubridate::days(16)

results <- NOAA_GEFS(
  outfolder,
  lat.in,
  lon.in,
  sitename,
  tz = NULL,
  start_date = start_date,
  end_date = end_date,
  overwrite = FALSE,
  verbose = FALSE
)

save(results, 
     file = paste0(outfolder, "/results.RData"))
cat("===== Cary Institute downloaded =====\n")

#------------------------
#    Harvard Forest 
#------------------------
lat.in <- 42.5369
lon.in <- -72.17266
sitename <- "HarvardForest"
outfolder <- "../GEFS/HarvardForest"
tz <- "America/New_York"
end_date <- as.POSIXct(start_date, tz = tz) + lubridate::days(16)

results <- NOAA_GEFS(
  outfolder,
  lat.in,
  lon.in,
  sitename,
  tz = NULL,
  start_date = start_date,
  end_date = end_date,
  overwrite = FALSE,
  verbose = FALSE
)
cat("===== Harvard Forest downloaded =====\n")



cat("\n ----- END ----- \n\n")
