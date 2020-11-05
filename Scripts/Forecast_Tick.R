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

use.gefs <- TRUE # do we use gefs (if false we climate ensembles)
weather.dir <- ifelse(use.gefs, "GEFS", "Climate")
latency <- 31    # data latency

library(rjags)
library(tidyverse)
library(lubridate)
library(ncdf4)
library(plantecophys)
library(rvest)
library(xlsx)

source("Functions/gefs_prior_functions.R")
source("Functions/ticks_2020.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/convergence_check.R")
source("Functions/scale_met_forecast.R")
source("Functions/cary_weather_ensembles.R")
source("Functions/create_ncdf_tick.R")
source("Models/tickForecastFilter_monthEffectLifeStage.R")

date.pattern <- "\\d{4}-\\d{2}-\\d{2}"

# =================================================== #
#                Site and dir set-up                  #
# =================================================== #

site.vec <- c("Green Control", "Henry Control", "Tea Control")

## read array job number and subset site
array.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# array.num <- 2 # for testing 
site.name <- site.vec[array.num]

cat("Running forecast on", site.name, "\n")

## parameter estimates from historical fit
top.dir <- "../FinalOut/A_Correct"   # all fitted models start here
model.type <- "RhoModels/RhoAllStart/MonthEffect/DiffPriors" # path to specific model 
site.dir.name <- gsub(" Control", "", site.name) # site directory name
rdata.name <- paste0("Combined_thinMat_RhoStartMonthEffect_LifeStageDiffPriors_", 
                     gsub(" ", "", site.vec[array.num]), 
                     ".RData") 
load(file.path(top.dir,
               model.type,
               site.dir.name,
               rdata.name))

# forecast output directory
specific.name <- "LifeStageMonthEffect_UpdateAllParameters"
latency.dir <- paste0("dataLatency_", latency)
out.dir <- file.path("../FinalOut/Forecast2020Cary", 
                     site.dir.name,
                     model.type,
                     specific.name,
                     weather.dir,
                     latency.dir) 

# parameter estimates in data list, only Apr:Dec fitted during training
data <- update_data(params.mat, 4:12) 

# 2020 observations
ticks <- ticks_2020(site.name)
obs.dates <- pull(ticks, Date)
obs <- ticks %>% 
  select(c(Larvae, Nymphs, Adult.total)) %>% 
  t()

# data latency, days after observation data is "posted"
posted.days <- obs.dates + latency

cat("Observations on:", as.character(obs.dates), "\n")
cat("Posted on:", as.character(posted.days), "\n")

# indexing for known met and gefs
n.adapt <- NULL
n.chains <- 5
n.iter <- 50000

forecast.start.day <- "2020-05-19" # day forecast stems from
every.day <- seq.Date(ymd(forecast.start.day), today(), by = 1) # everyday

# =================================================== #
#                Initial Conditions                   #
# =================================================== #
# need to set initial conditions for the first forecast
# will use the mean number of each life stage observed in 
# May from 1995-2018. Those observations are split among 
# two data sets
# ticks 1995-2005
train.ticks <- read.csv("/projectnb/dietzelab/fosterj/Data/tick_cleaned") %>% 
  filter(Grid == site.name) %>% 
  mutate(month = month(DATE)) %>% 
  filter(month == 5) %>% 
  select(c(n_larvae, n_nymphs, n_adults))

# ticks 2006-2018
hind.ticks <- get_ticks_2006_2018(site.name)
ticks.06.18 <- t(hind.ticks$obs) %>% 
  as.data.frame() %>% 
  mutate(DATE = (hind.ticks$date)) %>% 
  mutate(month = month(DATE)) %>% 
  filter(month == 5) %>% 
  select(c(n_larvae, n_nymphs, n_adults))

# bind together and calculate mean
ic.ticks <- bind_rows(train.ticks, ticks.06.18) 
ic.ticks <- round(apply(ic.ticks, 2, mean))

data$ic <- ic.ticks # initial condition, put in data for jags
y.vec <- rep(NA, 3) # no Ys at beginning of forecast, assimilated later

# =================================================== #
#       Observed met from Cary, previous months       #
# =================================================== #

met.cary <- read.csv("/projectnb/dietzelab/fosterj/Data/CaryCurrentYearMet.csv")
met.cary$DATE <- mdy(met.cary$DATE)
met.cary <- met.cary %>% 
  select(c("DATE","MAX_TEMP", "MIN_TEMP", "MAX_RH", "MIN_RH", "tot_prec")) 
end.date.met.cary <- last(met.cary$DATE)
 
# =================================================== #
#       Observed met from Cary, current month         #
# =================================================== #

# scrape current month table from Cary's website
month.url <- "http://www.caryinstitute.org/science/research-projects/environmental-monitoring-program/weather-climate"
webpage <- read_html(month.url)
tbls <- html_nodes(webpage, "table")

current.month <- webpage %>%
  html_nodes("table") %>%
  .[1] %>%                                # first table on the page is daily summaries
  html_table(fill = TRUE) %>% 
  as.data.frame() %>% 
  select(c("DATE", "MAX..TEMP...F.", "MIN..TEMP...F.", "TOTAL.PRECIP...in..")) %>% 
  rename("MAX_TEMP" = "MAX..TEMP...F.",   # rename to match met.cary
         "MIN_TEMP" = "MIN..TEMP...F.", 
         "tot_prec" = "TOTAL.PRECIP...in..") %>% 
  mutate("DATE" = mdy(DATE),
         "MAX_TEMP" = (MAX_TEMP-32)*5/9,  # F to C
         "MIN_TEMP" = (MIN_TEMP-32)*5/9,
         "MAX_RH" = NA,                   # don't have RH for current month
         "MIN_RH" = NA,
         "tot_prec" = tot_prec/25.4)  %>%     # in to mm
  filter(DATE > end.date.met.cary) # make sure we remove duplicate observations

# combine to single data set
met.cary <- bind_rows(met.cary, current.month)

gdd <- rep(NA, nrow(met.cary))
for(t in 1:nrow(met.cary)){
  gdd[t] <- max(mean(met.cary$MAX_TEMP[t], met.cary$MIN_TEMP[t]) - 10, 0)
}

# create vpd and cumulative gdd
met.cary$vpd <- RHtoVPD(met.cary$MIN_RH, met.cary$MIN_TEMP)
met.cary$cdd <- cumsum(gdd)

# scale to the training data means
hist.means <- scale_met_forecast()     # training met data means
met.cary$MAX_TEMP <- met.cary$MAX_TEMP - hist.means["MAX_TEMP"]
met.cary$MIN_TEMP <- met.cary$MIN_TEMP - hist.means["MIN_TEMP"]
met.cary$MAX_RH <- met.cary$MAX_RH - hist.means["MAX_RH"]
met.cary$MIN_RH <- met.cary$MIN_RH - hist.means["MIN_RH"]
met.cary$vpd <- met.cary$vpd - hist.means["vpd"]


# =================================================== #
#                   NOAA GEFS                         #
# =================================================== #
# GEFS has come in four different ways, using rnoaa,
# noaaGEFSpoint, noaaGEFSgrid, and from the archive system
# all have different directories, days associated

# directories
dir.gefs.rnoaa <- "/projectnb/dietzelab/fosterj/Data/GEFSrnoaa/Cary/" # from rnoaa downloads
dir.gefs.point <- "/projectnb/dietzelab/fosterj/Data/GEFSpoint/NOAAGEFS_6hr/CARY/" # from noaaGEFSpoint
dir.gefs.grid <- "/projectnb/dietzelab/fosterj/Data/GEFSgrid/NOAAGEFS_6hr/CARY/" # from noaaGEFSgrid
dir.gefs.arch <- "/projectnb/dietzelab/fosterj/Data/GEFSarchive/" # archived GEFS downloads

# get dates from the different directories

# rnoaa dates
files.rnoaa <- list.files(dir.gefs.rnoaa)
dates.rnoaa <- str_extract(files.rnoaa, date.pattern) %>% 
  ymd() %>% 
  discard(is.na)

# gefs point dates
files.point <- list.files(dir.gefs.point)
dates.point <- ymd(files.point)

# gefs grid dates
files.grid <- list.files(dir.gefs.grid)
dates.grid <- ymd(files.grid) 

# gefs archive dates
files.arch <- list.files(dir.gefs.arch, pattern = ".csv")
dates.arch <- str_extract(files.arch, "^\\d{8}") %>% 
  ymd()




# =================================================== #
#                    FORECAST                         #
# =================================================== #



for(i in seq_along(every.day)){
# i=1
  cat("\n\n===================================================\n\n")
  
  per <- i / length(every.day) * 100
  cat(round(per), "% completed\n")
  
  current.day <- every.day[i]  # date forecast is made
  
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
    
    # grab forecast from yesterday
    yesterday <- current.day - 1
    yest.dir <- file.path(out.dir, yesterday) # where forecasts are stored
    files.forecast <- list.files(yest.dir) 
    
    # assimilate forecast from climate ensembles
    if(any(grepl("clim", files.forecast))){
      files.forecast <- files.forecast[grep("clim", files.forecast)]
      nc <- nc_open(file.path(yest.dir, files.forecast[1])) # read nc file
      nrow <- nrow(ncvar_get(nc, "larvae")) # for dimensions
      params.nc <- matrix(0, 1, length(params.nc.names))
      colnames(params.nc) <- ncvar_get(nc, "parameter_names") # don't change
      nc.larva <- nc.nymph <- nc.adult <- matrix(NA, nrow, length(files.forecast))
      for(cf in seq_along(files.forecast)){
        nc <- nc_open(file.path(yest.dir, files.forecast[cf])) # read nc file
        f.days <- ncvar_get(nc, "time")     # dates for each forecast
        f.2.grab <- which(current.day == f.days) # the day we need
        nc.larva[,cf] <- ncvar_get(nc, "larvae")[,f.2.grab] # larva forecasts
        nc.nymph[,cf] <- ncvar_get(nc, "nymph")[,f.2.grab]  # nymph forecasts
        nc.adult[,cf] <- ncvar_get(nc, "adult")[,f.2.grab]  # adult forecasts
        params.ens <- ncvar_get(nc, "parameter_samples")
        colnames(params.ens) <- ncvar_get(nc, "parameter_names") # don't change
        params.nc <- rbind(params.nc, params.ens)
      }
      params.nc <- params.nc[-1,] # remove first row from initialization above
    } else {
      nc <- nc_open(file.path(yest.dir, files.forecast[1])) # read nc file
      f.days <- ncvar_get(nc, "time")     # dates for each forecast
      f.2.grab <- which(current.day == f.days) # the day we need
      nc.larva <- ncvar_get(nc, "larvae")[,f.2.grab] # larva forecasts
      nc.nymph <- ncvar_get(nc, "nymph")[,f.2.grab]  # nymph forecasts
      nc.adult <- ncvar_get(nc, "adult")[,f.2.grab]  # adult forecasts
      params.nc <- ncvar_get(nc, "parameter_samples")
      colnames(params.nc) <- ncvar_get(nc, "parameter_names") # don't change
    }
    
    # pull out median predictions for DA
    larva.m <- median(nc.larva)
    nymph.m <- median(nc.nymph)
    adult.m <- median(nc.adult)

    # update posteriors
    data <- update_data(params.nc, 1:12)	
    
    # initial condition is median of forecast
    data$ic <- c(larva.m, nymph.m, adult.m)
    
    # set forecast.start.day
    forecast.start.day <- new.start.day
    
  } 

  # subset observed met
  build.obs.met <- met.cary %>% 
    filter(DATE <= current.day) %>% 
    filter(DATE >= forecast.start.day)
  
  # cumulative growing degree days at end of observation period
  # needed for cum.gdd calculations across GEFS 
  end.cum.gdd <- build.obs.met$cdd[nrow(build.obs.met)]
  
  # grab appropriate GEFS
  gefs.needed <- ymd(current.day)+1 # need gefs starting 'tomorrow'
  if(use.gefs){
    use.climate <- FALSE
    cat("Weather forecast start day:", as.character(gefs.needed), "\n")
    if(gefs.needed %in% dates.grid){
      cat("Using GEFS from grid\n")
      met.gefs <- get_gefs_noaaGEFSpoint("grid", gefs.needed, end.cum.gdd)
    } else if(gefs.needed %in% dates.point){
      cat("Using GEFS from point\n")
      met.gefs <- get_gefs_noaaGEFSpoint("point", gefs.needed, end.cum.gdd)
    } else if(gefs.needed %in% dates.arch){
      cat("Using GEFS from archives\n")
      met.gefs <- get_gefs_archive(gefs.needed, end.cum.gdd)
    } else if(gefs.needed %in% dates.rnoaa){
      cat("Using GEFS from rnoaa\n")
      met.gefs <- get_gefs_rnoaa(gefs.needed, end.cum.gdd)
    } else {
      cat("No GEFS - using climate\n")
      use.climate <- TRUE
      clim.ens <- cary_weather_ensembles(yday(gefs.needed))
    }  
  } else {
    cat("Using climatology ensembles\n")
    use.climate <- TRUE
    clim.ens <- cary_weather_ensembles(yday(gefs.needed))
    # met.gefs <- met.gefs[sample(length(met.gefs), 7)]
  }
  
  if(use.climate){
    runs <- length(clim.ens) # the number of times we need to run the forecast - one per ensemble
    data$cum.gdd <- pull(build.obs.met, cdd)
 
    # total days in forecast
    data$n.days <- nrow(build.obs.met) + nrow(clim.ens[[1]])
    all.days <- seq.Date(ymd(forecast.start.day), length.out = data$n.days, by = 1)
  } else {
    
    runs <- 1 # the number of times we need to run the forecast - once
    
    # extract the mean vector and precision matrix for each weather variable
    cum.gdd <- get_gefs_mean_prec(met.gefs, "cum.gdd") 
    min.temp <- get_gefs_mean_prec(met.gefs, "min.temp", hist.means["MIN_TEMP"])
    
    # total days in forecast
    data$n.days <- nrow(build.obs.met) + length(min.temp$mu)
    all.days <- seq.Date(ymd(forecast.start.day), length.out = data$n.days, by = 1)
    
    # control flow for cumulative gdd priors if we drop 1st cum.gdd forecast
    if(!is.null(cum.gdd$add.2.obs)){
      data$cum.gdd <- c(pull(build.obs.met, cdd), cum.gdd$add.2.obs)
      
      # need day of year for gefs prior
      doys <- yday(all.days[-(1:(nrow(build.obs.met)+1))])
      gefs.priors <- gefs_prior(doys)
      data$cum.gdd.prior <- pull(gefs.priors, cum.gdd)
    } else {
      data$cum.gdd <- pull(build.obs.met, cdd)  
      
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
    
    
    data$cum.gdd.gefs.mu <- cum.gdd$mu
    data$cum.gdd.gefs.prec <- cum.gdd$prec
    data$obs.temp.gefs.mu <- min.temp$mu
    data$obs.temp.gefs.prec <- min.temp$prec
  }
  
  # build data needed for jags
  data$obs.temp <- pull(build.obs.met, MIN_TEMP)
  data$month.index <- month(all.days)
  
  # build y matrix
  y <- matrix(NA, 3, data$n.days)
  y[,1] <- y.vec # does not change unless we assimilate
  data$y <- y
  
  cat("Forecast initial conditions and start day:", as.character(forecast.start.day), "\n")
  cat("Forecast is being run on:", as.character(current.day), "\n")
  
  for(r in 1:runs){
    
    # control flow for file name and jags data for climate ensembles
    if(use.climate){
      cat("=== Fitting climate ensemble", r, "===\n")
      data$clim.gdd.mu <- clim.ens[[r]] %>% pull(cum.gdd) 
      data$clim.met.obs.mu <- clim.ens[[r]] %>% pull(min.temp) - hist.means["MIN_TEMP"]
      
      outname <- paste(forecast.start.day, 
                       current.day, 
                       paste0("clim", r),
                       "Tick_forecast.nc", 
                       sep = "_") 
    } else {
      outname <- paste(forecast.start.day, 
                       current.day, 
                       "Tick_forecast.nc", 
                       sep = "_") 
    }
    
    # compile and sample
    mcmc.model <- run_jagsFilter(data, 
                                 n.adapt, 
                                 n.chains, 
                                 n.iter, 
                                 use.climate)
    
    # check convergence
    out <- convergence_check(jags.out = mcmc.model$jags.out,
                             model = mcmc.model$compiled.model,
                             monitor = mcmc.model$monitor,
                             n.iter = n.iter/2,
                             min.eff.size = 2000)
    
    ## split output
    # jags.out <- mcmc.model$jags.out
    # out <- list(params = NULL, predict = NULL)
    # mfit <- as.matrix(jags.out, chains = TRUE)
    # pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
    # chain.col <- which(colnames(mfit) == "CHAIN")
    # out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
    # out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
    
    # extract parameters and predicted states
    params <- as.matrix(out$params)
    preds <- as.matrix(out$predict)
    
    # create and check dir, each day the forecast is initiated gets own dir
    out.dir.current.day <- file.path(out.dir, current.day)
    if(!dir.exists(out.dir.current.day)) dir.create(out.dir.current.day, recursive = TRUE)
    
    # full path
    ncfname <- file.path(out.dir.current.day, outname)
    
    # vector for what days tick data is assimilated (will need updating)
    data.assimilation <- as.numeric(complete.cases(t(data$y)))
    
    # write netcdf 
    create_ncdf_tick(ncfname = ncfname,
                     preds = preds,
                     params = params,
                     start.date = forecast.start.day,
                     data.assimilation = data.assimilation)

  }
  
  
  
  ## save forecasts within specific folders for each forecast date? or same directory?
}




