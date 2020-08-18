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

source("Functions/gefs_prior_functions.R")

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
for(i in seq_along(files)){
  met.ncdf <- nc_open(paste0(dir.obs.met, files[i]))  
  for(v in seq_along(var.names)){
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])
    met.observed[i, v] <- ncvar_get(met.ncdf, var.names[v])  
  }
}
colnames(met.observed) <- var.names

# last date of weather observations
last.met.obs.date <- ymd(met.observed.dates[length(met.observed.dates)])

# calculate growing degree days with base 10
base <- 10
gdd <- rep(NA, nrow(met.observed))
for(i in 1:nrow(met.observed)){
  gdd[i] <- max(mean(met.observed[i, "TMAX"], met.observed[i, "TMIN"]) - base, 0)
}

# add cumulative growing degree days to observed met
met.observed <- met.observed %>% 
  as.data.frame() %>% 
  mutate(cum.gdd = cumsum(gdd))

# cumulative growing degree days at end of observation period
# needed for cum.gdd calculations across GEFS 
end.cum.gdd <- met.observed$cum.gdd[nrow(met.observed)]

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

met.gefs <- list()
for(ens in 1:21){ # get data for each ensemble member
  gefs.ens <- nc_open(paste0(dir.gefs, days.2.grab, "/", ens.files[ens])) 
  gefs.data <- matrix(NA, n.days, n.var)
  for(v in seq_along(var.names)){
    gefs.data[,v] <- ncvar_get(gefs.ens, var.names[v])  
  }
  colnames(gefs.data) <- var.names
  
  gdd <- rep(NA, nrow(gefs.data))
  for(i in 1:nrow(gefs.data)){
    gdd[i] <- max(mean(gefs.data[i, "max.temp"], gefs.data[i, "min.temp"]) - base, 0)
  }
  
  gefs.cum.gdd <- cumsum(c(end.cum.gdd, gdd))
  gefs.data <- gefs.data %>% 
    as.data.frame() %>% 
    mutate(cum.gdd = gefs.cum.gdd[-1])
  
  met.gefs[[ens]] <- gefs.data
}

min.temp <- get_gefs_mean_prec(met.gefs, "min.temp")
max.temp <- get_gefs_mean_prec(met.gefs, "max.temp")
min.rh <- get_gefs_mean_prec(met.gefs, "min.rh")
max.rh <- get_gefs_mean_prec(met.gefs, "max.rh")
precip <- get_gefs_mean_prec(met.gefs, "precip")
vpd <- get_gefs_mean_prec(met.gefs, "vpd")
cum.gdd <- get_gefs_mean_prec(met.gefs, "cum.gdd")

# =================================================== #
#                     Testing                         #
# =================================================== #

## parameter estimates from last forecast
load("../FinalOut/A_Correct/RhoModels/RhoAllStart/MonthEffect/Henry/Combined_thinMat_RhoStartMonthEffectLifeStage_HenryControl.RData")

data <- update_data(params.mat, 4:12)

## initial conditions
data$ic <- rpois(3, c(500, 50, 10))
 
ticks <- get_ticks_2006_2018(grid)


# dates
ic.date <- last.met.obs.date - 15

met.obs.index <- match(c(ic.date+1, last.met.obs.date), ymd(met.observed.dates))
met.subset <- met.observed[met.obs.index[1]:met.obs.index[2],]

data$cum.gdd <- pull(met.subset, cum.gdd)
data$obs.temp <- pull(met.subset, TMIN)
data$n.days <- nrow(met.subset) + length(cum.gdd$mu)
data$seq.days <- (data$n.days - 1):1
data$cum.gdd.gefs.mu <- cum.gdd$mu
data$cum.gdd.gefs.prec <- cum.gdd$prec
data$obs.temp.gefs.mu <- min.temp$mu
data$obs.temp.gefs.prec <- min.temp$prec
data$n.days.gefs <- length(cum.gdd$mu)

# need the month of all days in the sequence
all.days <- seq.Date(ic.date, length.out = data$n.days, by = 1)
data$month.index <- month(all.days)

# need day of year for gefs
doys <- yday(all.days[-(1:nrow(met.subset))])
gefs.priors <- gefs_prior(doys)
data$cum.gdd.prior <- pull(gefs.priors, cum.gdd)
data$obs.temp.prior <- pull(gefs.priors, min.temp)


# indexing for known met and gefs
n.adapt <- 20
n.chains <- 3

j.model <- jags.model(
  file = textConnection(model),
  data = data,
  inits = inits,
  n.adapt = n.adapt,
  n.chains = n.chains
)

jags.out <- coda.samples(model = j.model,
                         variable.names = monitor,
                         n.iter = 50000)

jags.mat <- as.matrix(jags.out)
met.obs.jags <- jags.mat[,grep("met.obs.mu", colnames(jags.mat), fixed = TRUE)]
met.obs.jags.ci <- apply(met.obs.jags, 2, quantile, c(0.025, 0.05, 0.975))

par(mfrow=c(1,1))
plot(dat[1,], type = "l", ylim = c(0,20))
for(m in 2:nrow(dat)) lines(dat[m,])
for(m in 1:nrow(met.obs.jags.ci)) lines(met.obs.jags.ci[m,], col = "blue")

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

# need observed day indicies and gefs day indicies

