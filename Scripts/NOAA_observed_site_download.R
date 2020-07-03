# need to set this for CRON, will default to geo home instead of projectnb
setwd("/projectnb/dietzelab/fosterj/NEFI_tick")

library(rnoaa)
library(ncdf4)
library(dplyr)
library(lubridate)


source("Functions/NOAA_observed.R")

station.id <- "GHCND:USW00064756" # Millbrook, NY
sitename <- "CaryInstitute"

# site level storage
outfolder <-  "../NOAA_Stations/Cary"

cat("===== Attempting", sitename, "download =====\n")
NOAA_observed(outfolder, sitename, station.id)
cat("=====", sitename, "downloaded =====\n")




# for getting more than one day
# day.vec <- seq.Date(ymd("2020-04-01"), ymd("2020-06-26"), by = 1)

# for(i in seq_along(day.vec)){
#   NOAA_observed(outfolder, sitename, station.id, day.2.grab = day.vec[i])
# }




dat <- nc_open("../NOAA_Stations/Cary/NOAA_Observed.CaryInstitute.2020-06-26.nc")
ncvar_get(dat, "TMIN")




