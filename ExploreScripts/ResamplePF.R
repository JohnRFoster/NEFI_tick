library(ecoforecastR)

## Load Model Function
source("Models/ensemble_gdd_k_window.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/get_last_day.R")
source("Functions/Future_Met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/life_stage_ci.R")
source("Functions/site_data_met.R")
dir <- "../FinalOut/HB_Partial_GDD"

## Specifications
model.type <- "K_estimate"
version <- "WindowLoop_4alpha"
Rdata <- "GDDSwitch_K_Window4alpha_AllChains.RData"

site <- "Green Control"
if(site == "Green Control"){
  s <- 1
} else if(site == "Henry Control"){
  s <- 2
} else {
  s <- 3
}

# number of ensembles to run
Nmc <- 10

## Load Model Fit Output
load(file.path(dir, model.type, version, Rdata))
params <- as.matrix(out.test$params)
draw <- sample.int(nrow(params), Nmc, replace = TRUE)
params <- params[draw,]
params.orig <- params

## probability of detection
theta.larvae <- params[, "theta.larvae"]
theta.nymph <- params[, "theta.nymph"]
theta.adult <- params[, "theta.adult"]

## Load Future Met
met <- future_met(site)

## Load Future ticks
ticks.future <- get_ticks_2006_2018(site)

ensemble.out <- pf.out <- list()

N_days <- 16 # number of days to forecast

# load historical data and reduce to site
data <- cary_ticks_met_JAGS()
hist.data <- site_data_met(site, NULL, data)

# first initial condition is last obs in historical timeseries
last.obs <- hist.data$y[,ncol(hist.data$y)]
ic <- matrix(NA, 3, Nmc)
for(i in 1:3){
  # create initial condition ensemble members
  ic[i,] <- rpois(Nmc, last.obs) 
}

## all days in future sequence
last.day <- get_last_day()
last.day <- as.Date(last.day[1,2])
future.days <- seq.Date(last.day + 1, last.day + (365*8), by = 1)

## index of when observations happen (days from last date)
obs.index <- vector()
vec <- which(ticks.future$date %in% future.days)
for(i in seq_along(which(ticks.future$date %in% future.days))){
  xx <- which(ticks.future$date[i] == future.days)
  obs.index <- c(obs.index, xx)
}

day <- 0
resamp.archive <- list()
output <- array(NA, dim = c(3, # number of forecasts
                            length(future.days), # life stages
                            Nmc)) # number of ensemble members

update_params <- function(params,index){
  params <- params[index,]
  return(params)
}

# for(t in seq_along(future.days)){
for(t in 1:320){
  ## month counter
  if(t %% 31 == 0){print(t)}
  
  ## run ensemble
  gdd <- met$cum.gdd[t:(N_days + t - 1)]
  out  <- ensemble_gdd_k_window(site, params, ic, gdd, N_days, Nmc)
  resamp.archive[[t]] <- out
  
  ## initial condition for next forecast
  ic[1,] <- out[1, 1,]
  ic[2,] <- out[2, 1,]
  ic[3,] <- out[3, 1,]
}
  ## analysis step - if observation is present
  if(t %in% obs.index){
    day <- day + 1
    
    cat("Observation", day, "\n")
    
    # reset N_days
    build <- N_days
    
    # first day forecast ensembles
    fore <- resamp.archive[[t-N_days]][,build,]
    
    # loop over days that have forecasts
    for(d in t+1-N_days:2){
      build <- build-1
      fore <- cbind(fore, resamp.archive[[d]][,build,])
    }
    
    like <- apply(fore, 2, mean)
    # like <- mean(fore)
    
    wt <- dpois(ticks.future$obs[, day], like)
    # hist(wt)
    
    ## resample ensembles 
    N.samp <- length()
    index <- sample.int(Nmc, Nmc, replace = TRUE, prob = wt)
    
    ## initial condition for next forecast
    ic <- output[, t, index]
    
    ## initial condition for next forecast if we observe 
    # ic[1,] <- rpois(Nmc, ticks.future$obs[,day])
    # ic[2,] <- rpois(Nmc, ticks.future$obs[,day])
    # ic[3,] <- rpois(Nmc, ticks.future$obs[,day])
    
    ## resample parameters
    params <- update_params(params, index)
    params.archive[[day-1]] <- params
  }
}

ReSampPF <- output
save(ReSampPF, file = "ReSampPF_Test_Green.RData")