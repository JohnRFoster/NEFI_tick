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
source("Functions/scale_met_forecast.R")
source("Functions/create_ncdf_tick.R")
source("Models/tickForecastFilter_monthEffectLifeStage.R")

date.pattern <- "\\d{4}-\\d{2}-\\d{2}"

# =================================================== #
#                  Testing set up                     #
# =================================================== #

site.vec <- c("Green Control", "Henry Control", "Tea Control")

## read array job number and subset site
# array.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
array.num <- 2 # for testing
site.name <- site.vec[array.num]

cat("Running forecast on", site.name, "\n")

## training met data means
hist.means <- scale_met_forecast()

## parameter estimates from historical fit
top.dir <- "../FinalOut/A_Correct"   # all fitted models start here
model.type <- "RhoModels/RhoAllStart/MonthEffect" # path to specific model 
site.dir.name <- gsub(" Control", "", site.name) # site directory name
rdata.name <- "Combined_thinMat_RhoStartMonthEffectLifeStage_HenryControl.RData" 
load(file.path(top.dir,
               model.type,
               site.dir.name,
               rdata.name))

# forecast output directory
specific.name <- "LifeStageMonthEffect_UpdateAllParameters"
today.dir <- paste0("Submit_", today())
out.dir <- file.path("../FinalOut/ForecastTestRuns", today.dir, "Cary", 
                     model.type,
                     specific.name,
                     site.dir.name)

# check and create dir
if(!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

# parameter estimates in data list, only Apr:Dec fitted during training
data <- update_data(params.mat, 4:12) 

# using 2012 time series to test work flow
ticks <- get_ticks_2006_2018("Henry Control")
index.2012 <- which(year(ticks$date) == 2012)
obs <- ticks$obs[,index.2012]   # 2012 observations
y.vec <- obs[,1]
obs <- obs[,2:5] # remove first set
data$ic <- rpois(3, y.vec+1)         # initial condition, put in data for jags

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
# need to do this before centering other variables to historical means
base <- 10
gdd <- rep(NA, nrow(met.observed))
for(i in 1:nrow(met.observed)){
  gdd[i] <- max(mean(met.observed[i, "TMAX"], met.observed[i, "TMIN"]) - base, 0)
}

# add cumulative growing degree days to observed met
met.observed <- met.observed %>% 
  as.data.frame() %>% 
  mutate(cum.gdd = cumsum(gdd)) %>% 
  mutate(Date = ymd(met.observed.dates)) %>% 
  mutate(TMAX = TMAX - hist.means["MAX_TEMP"]) %>% # center to historical means
  mutate(TMIN = TMIN - hist.means["MIN_TEMP"])

dir.gefs <- "../GEFS/Cary/"
files.gefs <- list.files(dir.gefs)

met.gefs.dates <- ymd(str_extract(files.gefs, date.pattern))

# indexing for known met and gefs
n.adapt <- 5000
n.chains <- 5
n.iter <- 100000

# forecast loop, need to check:
# 1. if there is a new tick observation
# 2. the latest observed weather
every.day <- seq.Date(met.gefs.dates[1]+4, today(), by = 1) # everyday 
forecast.start.day <- "2020-05-20" # day forecast stems from

for(i in seq_along(every.day)){
# i=1
  cat("\n\n===================================================\n\n")
  
  per <- i / length(every.day) * 100
  cat(round(per), "percent completed\n")
  
  current.day <- every.day[i]  # date forecast is made
  
  # usually run three or four days behind
  last.met.obs.date <- current.day - 4 
  
  # need to check if there is a new tick observation
  # this determines new initial condition and forecast start day
  if(current.day %in% posted.days){
    # match and assimilate
    index <- which(current.day == posted.days)
    y.vec <- obs[,index] # update 
    
    # new date forecast stems from
    new.start.day <- obs.dates[index]
    
    cat("New tick observation on", as.character(obs.dates[index]), "\n")
    cat("Posted on", as.character(posted.days[index]), "\n")
    
    # grab forecast; all files written
    files.forecast <- list.files(out.dir) 
    
    # the files that stem from previous forecast.start.day
    forecast.group <- grep(forecast.start.day, files.forecast) 
    files.forecast <- files.forecast[forecast.group] # subset
    assim.file <- files.forecast[length(files.forecast)]
    
    load(file.path(out.dir, assim.file))
    
    # update posteriors 
    data <- update_data(parameters, 1:12) 
    
    # pull out median predictions
    forecast.med <- round(apply(preds, 2, median))
    l.seq <- seq(1, by = 3, length.out = ncol(preds)/3) # larva columns
    
    larva <- forecast.med[l.seq]
    nymph <- forecast.med[l.seq+1] # add 1 for nymph columns
    adult <- forecast.med[l.seq+2] # add 2 for adult columns
    
    # find day that we need
    median.index <- which(date.fore == new.start.day)
    
    # initial condition is median of forecast
    data$ic <- c(larva[median.index], nymph[median.index], adult[median.index])
    
    # set forecast.start.day
    forecast.start.day <- new.start.day
    
  } 

  # subset observed met
  build.obs.met <- met.observed %>% 
    filter(Date <= last.met.obs.date) %>% 
    filter(Date >= forecast.start.day)
  
  # cumulative growing degree days at end of observation period
  # needed for cum.gdd calculations across GEFS 
  end.cum.gdd <- build.obs.met$cum.gdd[nrow(build.obs.met)]
  
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
  days.2.grab <- files.gefs[gefs.needed] # gefs date we need
  
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
  cum.gdd <- get_gefs_mean_prec(met.gefs, "cum.gdd") 
  min.temp <- get_gefs_mean_prec(met.gefs, "min.temp", hist.means["MIN_TEMP"])
  
  # total days in forecast
  data$n.days <- nrow(build.obs.met) + length(min.temp$mu)
  
  # need the month of all days in the sequence
  all.days <- seq.Date(ymd(forecast.start.day), length.out = data$n.days, by = 1)
  data$month.index <- month(all.days)
  
  # control flow for cumulative gdd priors if we drop 1st cum.gdd forecast
  if(!is.null(cum.gdd$add.2.obs)){
    data$cum.gdd <- c(pull(build.obs.met, cum.gdd), cum.gdd$add.2.obs)
    
    # need day of year for gefs prior
    doys <- yday(all.days[-(1:(nrow(build.obs.met)+1))])
    gefs.priors <- gefs_prior(doys)
    data$cum.gdd.prior <- pull(gefs.priors, cum.gdd)
  } else {
    data$cum.gdd <- pull(build.obs.met, cum.gdd)  
    
    # need day of year for gefs prior
    doys <- yday(all.days[-(1:nrow(build.obs.met))])
    gefs.priors <- gefs_prior(doys)
    data$cum.gdd.prior <- pull(gefs.priors, cum.gdd)
  }
  
  # need day of year for gefs
  doys <- yday(all.days[-(1:nrow(build.obs.met))])
  gefs.priors <- gefs_prior(doys)
  data$obs.temp.prior <- pull(gefs.priors, min.temp)
  
  # number of days in each ensemble variable forecast
  data$n.gefs.gdd <- length(cum.gdd$mu)
  data$n.gefs.min.temp <- length(min.temp$mu)
  
  # build rest of data needed for jags
  data$obs.temp <- pull(build.obs.met, TMIN)
  data$cum.gdd.gefs.mu <- cum.gdd$mu
  data$cum.gdd.gefs.prec <- cum.gdd$prec
  data$obs.temp.gefs.mu <- min.temp$mu
  data$obs.temp.gefs.prec <- min.temp$prec
  
  # build y matrix
  y <- matrix(NA, 3, data$n.days)
  y[,1] <- y.vec # does not change unless we assimilate
  data$y <- y
  
  cat("Forecast initial conditions and start day:", as.character(forecast.start.day), "\n")
  cat("Forecast is being run on:", as.character(current.day), "\n")
  cat("Most recent weather observation:", as.character(last.met.obs.date), "\n")
  cat("The GEFS file assimilated:", days.2.grab, "\n")
  
  # compile and sample
  mcmc.model <- run_jagsFilter(data, n.adapt, n.chains, n.iter)
  
  # check convergence
  out <- convergence_check(jags.out = mcmc.model$jags.out,
                           model = mcmc.model$compiled.model,
                           monitor = mcmc.model$monitor,
                           n.iter = n.iter/2,
                           min.eff.size = 2000)
  
  # extract parameters and predicted states
  params <- as.matrix(out$params)
  preds <- as.matrix(out$predict)
  
  # save as .RData
  outname <- paste(forecast.start.day, current.day, "Tick_forecast_jagsFilter.nc", sep = "_") # name file
  ncfname <- file.path(out.dir, outname)
  
  # vector for what days tick data is assimilated (will need updating)
  data.assimilation <- as.numeric(complete.cases(t(data$y)))
  
  # write netcdf 
  create_ncdf_tick(ncfname = ncfname,
                   preds = preds,
                   params = params,
                   start.date = forecast.start.day,
                   data.assimilation = data.assimilation)
  
  save(preds, params, ncfname,
       file = file.path(out.dir, "test.params.RData"))
  
  ## save forecasts within specific folders for each forecast date? or same directory?
}




