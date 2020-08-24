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

library(ecoforecastR)
library(tidyverse)
library(lubridate)
library(ncdf4)

source("Functions/gefs_prior_functions.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/convergence_check.R")
source("Models/tickForecastFilter_monthEffectLifeStage.R")

date.pattern <- "\\d{4}-\\d{2}-\\d{2}"

# =================================================== #
#                  Testing set up                     #
# =================================================== #

## parameter estimates from last forecast
top.dir <- "../FinalOut/A_Correct"   # all fitted models start here
model.type <- "RhoModels/RhoAllStart/MonthEffect" # path to specific model 
site.dir.name <- "Henry" # site
rdata.name <- "Combined_thinMat_RhoStartMonthEffectLifeStage_HenryControl.RData" 
load(file.path(top.dir,
               model.type,
               site.dir.name,
               rdata.name))

# output directory
specific.name <- "LifeStageMonthEffect_UpdateAllParameters"
out.dir <- file.path("../FinalOut/ForecastTestRuns/Cary", 
                     model.type,
                     specific.name,
                     site.dir.name)

# parameter estimates in data list, only Apr:Dec fitted during training
data <- update_data(params.mat, 4:12) 

# using 2012 time series to test work flow
ticks <- get_ticks_2006_2018("Henry Control")
index.2012 <- which(year(ticks$date) == 2012)
obs <- ticks$obs[,index.2012]   # 2012 observations
obs <- obs[,2:5]
data$ic <- rpois(3, 20)         # initial condition, put in data for jags

# 2012 observation dates, change year to 2020
obs.dates <- ymd(format(ticks$date[index.2012], "2020-%m-%d"))   
obs.dates <- obs.dates[2:5] # first and last dates outside range
ic.date <- obs.dates[1]   # date of observation
latency <- rpois(length(obs.dates), 31)   # data latency, days after observation data is "posted"
posted.days <- obs.dates + latency

# =================================================== #
#                    Observed met                     #
# =================================================== #

dir.obs.met <- "../NOAA_Stations/Cary/" # where met is stored
files <- list.files(dir.obs.met)

# need to get number of variables to store
met.ncdf <- nc_open(paste0(dir.obs.met, files[1]))
n.var <- length(met.ncdf$var) # number of variables
var.names <- names(met.ncdf$var) # variable names
met.observed.dates <- str_extract(files, date.pattern) # dates with  observed met

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
  mutate(cum.gdd = cumsum(gdd)) %>% 
  mutate(Date = ymd(met.observed.dates))

# cumulative growing degree days at end of observation period
# needed for cum.gdd calculations across GEFS 
end.cum.gdd <- met.observed$cum.gdd[nrow(met.observed)]

# =================================================== #
#                      NOAA GEFS                      #
# =================================================== #

dir.gefs <- "../GEFS/Cary/"
files <- list.files(dir.gefs)

forecast.start.day <- ic.date + 1

met.gefs.dates <- ymd(str_extract(files, date.pattern))

met.obs.index <- match(c(ic.date+1, last.met.obs.date), ymd(met.observed.dates))
met.subset <- met.observed[met.obs.index[1]:met.obs.index[2],]




# indexing for known met and gefs
n.adapt <- 20
n.chains <- 3
n.iter <- 50000


# jags.mat <- as.matrix(jags.out)
# met.obs.jags <- jags.mat[,grep("met.obs.mu", colnames(jags.mat), fixed = TRUE)]
# met.obs.jags.ci <- apply(met.obs.jags, 2, quantile, c(0.025, 0.05, 0.975))

# par(mfrow=c(1,1))
# plot(dat[1,], type = "l", ylim = c(0,20))
# for(m in 2:nrow(dat)) lines(dat[m,])
# for(m in 1:nrow(met.obs.jags.ci)) lines(met.obs.jags.ci[m,], col = "blue")


# forecast loop, need to check:
# 1. if there is a new tick observation
# 2. the latest observed weather
every.day <- seq.Date(met.gefs.dates[1]+4, today(), by = 1) # everyday 
forecast.start.day <- "2020-05-20" # day forecast stems from
for(t in seq_along(every.day)){
  current.day <- every.day[t]  # date forecast is made
  
  # usually run three or four days behind
  last.met.obs.date <- current.day - 4 
  
  # need to check if there is a new tick observation
  # this determines new initial condition and forecast start day
  if(current.day %in% posted.days){
    # match and assimilate
    index <- which(current.day == posted.days)
    data$ic <- obs[,index]
    forecast.start.day <- obs.dates[1]
    
    
  } 

  # grab GEFS
  gefs.needed <- match(ymd(last.met.obs.date)+1, met.gefs.dates)
  
  # some gefs days did not download, 
  # this finds the most recent day for non-existent gefs files
  gefs.diff <- 0 # reset 
  if(is.na(gefs.needed)){
    # time difference between last observed date and dates we have gefs forecasts
    lag.days <- difftime(ymd(last.met.obs.date)+1, met.gefs.dates)
    gefs.needed <- which(lag.days == min(lag.days[lag.days >= 0], na.rm = TRUE))
    gefs.diff <- as.numeric(difftime(ymd(last.met.obs.date)+1, met.gefs.dates[gefs.needed]))
  }
  days.2.grab <- files[gefs.needed] # gefs date we need
  
  # dir contains the gefs forecasts where each ensemble member is a separate file
  gefs.forecast.dir <- paste0(dir.gefs, days.2.grab)
  
  # read the dir and compile ensembles, also calculate cumulative gdd
  # returns a list of data frames, one for each ensemble
  met.gefs <- get_gefs_ens(gefs.forecast.dir, end.cum.gdd)
  
  # if there has been a lapse in gefs
  # need to remove the rows for which we have observations
  if(gefs.diff > 0){
    met.gefs <- modify(met.gefs, function(x) x[-c(1:gefs.diff),])
    cat("There were missing GEFS files\n")
  }
  
  # extract the mean vector and precision matrix for each weather variable
  min.temp <- get_gefs_mean_prec(met.gefs, "min.temp")
  cum.gdd <- get_gefs_mean_prec(met.gefs, "cum.gdd")
  
  # subset observed met
  build.obs.met <- met.observed %>% 
    filter(Date <= last.met.obs.date) %>% 
    filter(Date >= forecast.start.day)
  
  data$cum.gdd <- pull(build.obs.met, cum.gdd)
  data$obs.temp <- pull(build.obs.met, TMIN)
  data$n.days <- nrow(build.obs.met) + length(cum.gdd$mu)
  data$seq.days <- (data$n.days - 1):1
  data$cum.gdd.gefs.mu <- cum.gdd$mu
  data$cum.gdd.gefs.prec <- cum.gdd$prec
  data$obs.temp.gefs.mu <- min.temp$mu
  data$obs.temp.gefs.prec <- min.temp$prec
  data$n.days.gefs <- length(cum.gdd$mu)
  
  # need the month of all days in the sequence
  all.days <- seq.Date(ymd(forecast.start.day), length.out = data$n.days, by = 1)
  data$month.index <- month(all.days)
  
  # need day of year for gefs
  doys <- yday(all.days[-(1:nrow(build.obs.met))])
  gefs.priors <- gefs_prior(doys)
  data$cum.gdd.prior <- pull(gefs.priors, cum.gdd)
  data$obs.temp.prior <- pull(gefs.priors, min.temp)
  
  cat("Forecast initial conditions and start day:", as.character(forecast.start.day), "\n")
  cat("Forecast is being run on:", as.character(current.day), "\n")
  cat("Most recent weather observation:", as.character(last.met.obs.date), "\n")
  cat("The GEFS file assimilated:", days.2.grab, "\n")
  
  # compile and sample
  mcmc.model <- run_jagsFilter(data, n.adapt, n.chains, n.iter)
  
  # check convergence
  out <- convergence_check(mcmc.model$jags.out, mcmc.model$compiled.model)
  
  # update posteriors 
  parameters <- as.matrix(out$params)
  preds <- as.matrix(out$predict)
  data <- update_data(parameters, 1:12) 
  
  # thin for saving
  thin <- seq(1, nrow(parameters), length.out = 10000)
  preds <- preds[thin,]
  parameters <- parameters[thin,]
  
  # save as .RData
  outname <- paste(current.day, "Tick_forecast_jagsFilter.RData", sep = "_") # name file
  save(preds, parameters, data,
       file = file.path(out.dir, outname))
  
}




