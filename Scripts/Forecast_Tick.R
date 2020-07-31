# =================================================== #
#   This script is for forecasting at a single site   #
#   at the the Cary Institute of Ecosystem Studies    #
#                                                     #
#   The method for data assimilation is a jags filter #
#                                                     #
#   Jags model code is called from it's own script    #
#                                                     #
#   Here, we use observed met from the Millbrook, NY  #
#   NOAA met station, then switch to NOAA GEFS for    #
#   the forecast period.                              #   
# =================================================== #

library(rjags)
library(tidyverse)
library(lubridate)
library(ncdf4)


date.pattern <- "\\d{4}-\\d{2}-\\d{2}"

# =================================================== #
#                    Observed met                     #
# =================================================== #

dir.obs.met <- "../NOAA_Stations/Cary/" # where met is stored
files <- list.files(dir.obs.met)

# need to get number of variables to store
met.ncdf <- nc_open(paste0(dir.obs.met, files[1]))
n.var <- length(met.ncdf$var) # number of variables
var.names <- names(met.ncdf$var) # variable names
met.observed.dates <- str_extract(files, date.pattern) # dates

met.observed <- matrix(NA, length(files), n.var) # met data
met.observed.dates <- rep(NA, length(files))     # met dates (observed)
for(i in seq_along(files)){
  met.ncdf <- nc_open(paste0(dir.obs.met, files[i]))  
  for(v in seq_along(var.names)){
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])  
  }
}
colnames(met.observed) <- var.names

last.met.obs.date <- ymd(met.observed.dates[length(met.observed.dates)])


# =================================================== #
#                      NOAA GEFS                      #
# =================================================== #

dir.gefs <- "../GEFS/Cary/"
files <- list.files(dir.gefs)

forecast.start.day <- last.met.obs.date + 1

met.gefs.dates <- ymd(str_extract(files, date.pattern))
gefs.needed <- which(met.gefs.dates == forecast.start.day)
days.2.grab <- files[gefs.needed]

# dir contains the gefs forecasts where each ensemble member is a separate file
ens.files <- list.files(paste0(dir.gefs, days.2.grab))
ens.files <- ens.files[grepl(".nc", ens.files)] # want only .nc files

# get dimensions
gefs.ens <- nc_open(paste0(dir.gefs, days.2.grab, "/", ens.files[1])) 
n.var <- length(gefs.ens$var)
var.names <- names(gefs.ens$var)
n.days <- length(ncvar_get(gefs.ens, var.names[1]))

gefs.data <- matrix(NA, n.days, n.var)
met.gefs <- list()
for(ens in 1:21){ # get data for each ensemble member
  gefs.ens <- nc_open(paste0(dir.gefs, days.2.grab, "/", ens.files[ens])) 
  for(v in seq_along(var.names)){
    gefs.data[,v] <- ncvar_get(gefs.ens, var.names[v])  
  }
  colnames(gefs.data) <- var.names
  met.gefs[[ens]] <- gefs.data
}

# function to get the mean and precision from gefs ensembles
get_gefs_mean_prec <- function(met.gefs, met.var){
  dat <- t(map_dfc(met.gefs, function(x) x %>% as.data.frame() %>% select(met.var)))
  mu <- apply(dat, 2, mean)
  prec <- 1 / apply(dat, 2, var)  
  prec[is.infinite(prec)] <- 1e10 # if no variance (all ensembles the same), use very high prec
  return(list(mu = mu, prec = prec))
}

min.temp <- get_gefs_mean_prec(met.gefs, "min.temp")
max.temp <- get_gefs_mean_prec(met.gefs, "max.temp")
min.rh <- get_gefs_mean_prec(met.gefs, "min.rh")
max.rh <- get_gefs_mean_prec(met.gefs, "max.rh")
precip <- get_gefs_mean_prec(met.gefs, "precip")
vpd <- get_gefs_mean_prec(met.gefs, "vpd")


# =================================================== #
#                 Initial Conditions                  #
# =================================================== #

# grab initial conditions from last forecast



# =================================================== #
#     Check for and grab latest tick observations     #
# =================================================== #

# need to figure out how to do this part (google sheets??)



# =================================================== #
#               Build data list for jags              #
# =================================================== #



