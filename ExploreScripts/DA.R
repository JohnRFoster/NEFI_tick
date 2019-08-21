## Data Assimilation Code

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
Nmc <- 5000

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

N_days <- 16 # number of days to forecast

# load historical data and reduce to site
data <- cary_ticks_met_JAGS()
hist.data <- site_data_met(site, NULL, data)

# first initial condition is last obs in historical timeseries
last.obs <- hist.data$y[,ncol(hist.data$y)]
ic <- matrix(NA, 3, Nmc)
for(i in 1:3){
  # create initial condition ensemble members
  ic[i,] <- rpois(Nmc, last.obs[i]) 
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

# obs.index <- c(1, obs.index)
day <- 0
pf.archive <- list()
tick.pf <- array(NA, dim = c(3, # number of forecasts
                          3, # life stages
                          length(obs.index))) # quantile = c(0.025, 0.5, 0.975)

for(t in seq_along(future.days)){
# for(t in 1:320){
# for(t in 1:663){
  ## month counter
  if(t %% 31 == 0){print(t)}
  
  ## run ensemble
  gdd <- met$cum.gdd[t:(N_days + t - 1)]
  out <- ensemble_gdd_k_window(site, params, ic, gdd, N_days, Nmc)
  pf.archive[[t]] <- out
  
  ## update ic
  ic[1,] <- out[1, 1,]
  ic[2,] <- out[2, 1,]
  ic[3,] <- out[3, 1,]

  ## if new observation
  if(t %in% obs.index){
    # observation counter
    day <- day + 1
    cat("Observation", day, "\n")
    
    # reset N_days
    build <- N_days
    
    # first day forecast ensembles
    fore <- pf.archive[[t-N_days]][,build,]
    
    # loop over days that have forecasts
    for(d in t+1-N_days:2){
      build <- build-1
      fore <- cbind(fore, pf.archive[[d]][,build,])
    }
    
    like.all <- matrix(NA, 3, ncol(fore))
    for(i in 1:ncol(fore)){ # loop over ensembles
      
      ## calculate log likelihoods
      tick.like <- dpois(ticks.future$obs[,day], fore[,i], log = TRUE) 
      # tick.like <- dpois(round(fore[,i]), ticks.future$obs[,day-1], log = TRUE) 
      
      ## missing data as weight 1; log(1) = 0
      tick.like[is.na(tick.like)] <- 0
      
      ## convert to cumulative likelihood
      like.all[,i] <- (exp(cumsum(tick.like)))
      # tick.like <- exp(cumsum(tick.like))
    }
    
    ## mean weight at each time point
    #wbar <- mean(like.all)  
    wbar <- apply(like.all, 1, mean)  

    ## calculate weighted median and CI
    for(ls in 1:3){
      tick.pf[ls,,day] <- wtd.quantile(fore[ls,], 
                              like.all[ls,] / wbar[ls],
                              c(0.025, 0.5, 0.975)) 
    }
    
    ## create ic
    ## binary outcome of median pf by life stage and
    ## blend the poisson and zero inflation models
    inflate.larvae <- tick.pf[1, 2, day]*rbinom(Nmc, 1, theta.larvae) + 0.1
    inflate.nymph <- tick.pf[2, 2, day]*rbinom(Nmc, 1, theta.nymph) + 0.1
    inflate.adult <- tick.pf[3, 2, day]*rbinom(Nmc, 1, theta.adult) + 0.1
    
    ## initial condition for next forecast if we observe 
    ic[1,] <- rpois(Nmc, inflate.larvae)
    ic[2,] <- rpois(Nmc, inflate.nymph)
    ic[3,] <- rpois(Nmc, inflate.adult)
  }
}

save(tick.pf, file = "nonReSampPF_Test_Green.RData")
save(pf.archive, file = "nonReSampPFarchive_Test_Green.RData")











