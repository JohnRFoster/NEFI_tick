library(LaplacesDemon)
library(plantecophys)
library(ecoforecastR)

script.start <- Sys.time()

source("Functions/site_data_met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/ua_parts.R")
source("Functions/obs_prob.R")
source("Functions/Future_Met.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/tick_forecast.R")

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
Forecast_id <- "601144a5-ed57-4166-b943-38bfd66f9081"
ForecastProject_id <- 20200519 # Some ID that applies to a set of forecasts

grid <- "Green Control" # full site name
grid.short <- gsub(" Control", "", grid) # site sub directory name
grid.no.space <- gsub(" ", "", grid) # site name without spaces

met.driver <- "vpd"

# load dir
in.dir <- "../FinalOut/A_Correct"
model.structure <- "ObsProcModels"
model.drivers <- "Obs_L1N0A0.Proc_AdultVPD" 
model.name <- paste0("Combined_thinMat_Obs_L1.N0.A0_Proc_AdultVPD_",
                     grid.no.space,
                     ".RData")

# load training fit
load(file.path(in.dir, 
               model.structure, 
               model.drivers, 
               grid.short, 
               model.name))

Nmc <- 100
draw <- sample.int(nrow(params.mat), Nmc, replace = TRUE)
params <- params.mat[draw,]
states <- predict.mat[draw,]
params.hist <- list() # for storing through time
params.hist[[1]] <- params # store original draws

# output directory
out.dir <- file.path("../FinalOut/DA_Runs/Tick_Hindcast", 
                     model.structure,
                     model.drivers,
                     grid.short)

param.name <- paste0(grid.short, model.drivers, "_parameters.RData")

## Load Future ticks
ticks <- get_ticks_2006_2018(grid)
days <- ticks$date # date of observations
days <- days[days <= "2017-12-31"] # run hindcast through 2017 only, 2018 met data is NA
ticks.observed <- ticks$obs[,1:length(days)] # observed ticks
obs.delts <- ticks$df # number of days between sampling events

future.obs.index <- cumsum(as.numeric(diff(days)))
future.obs.index <- c(1, future.obs.index + 1)

# load and subset future met
met <- future_met(grid)
met$date <- as.character(met$date)
start.date <- as.character(days[1])
start.index <- match(start.date, met$date)
met <- met[start.index:nrow(met), ]
met$vpd <- RHtoVPD(met$min.rh, met$min.temp)

# first initial conditions is first observation of 2006
ua <- ua_parts(params, states, c("parameter", "ic", "process")) # uncertainty partitioning
obs.prob <- obs_prob(ua, 1, met$min.temp.scale[1]) # observation probability for zero inflation
ic <- matrix(NA, Nmc, 3) # store
first.observation <- pmax(1, ticks.observed[,1]) # turn 0s to 1s for rpois draws

# zero inflated Poisson initial condition
ic[,1] <- rpois(Nmc, first.observation[1]) * rbinom(Nmc, 1, obs.prob$theta.larva)
ic[,2] <- rpois(Nmc, first.observation[2]) * rbinom(Nmc, 1, obs.prob$theta.nymph)
ic[,3] <- rpois(Nmc, first.observation[3]) * rbinom(Nmc, 1, obs.prob$theta.adult)

cat("Larva initial condition range:", range(ic[,1]), "\n")
cat("Nymph initial condition range:", range(ic[,2]), "\n")
cat("Adult initial condition range:", range(ic[,3]), "\n")

t <- 1 # for testing
like <- matrix(1, Nmc, length(days))
resample <- rep(NA, length(days))
weight.ncdf <- rep(1, Nmc) # weights for ncdf
index <- 1
check.day <- rep(NA, length(days))

# for (t in 1:(length(days)-1)) {
for (t in 1:10) {
  
  ### forecast step ###
  forecast_issue_time <- as.Date(days[t])
  forecast.start.day <- as.character(days[t])  # date forecast issued
  forecast.end.day <- as.character(days[t+1])  # next observation date
  check.day[t] <- as.integer(days[t+1] - days[t]) # day we evaluate forecast and run DA
  check.day[t] <- max(2, check.day[t]) # one-day forecasts evaluate 2 index as first is ic
  
  # always make at least 16-day forecast
  if(check.day[t] < 16)  forecast.end.day <- as.character(days[t] + 16)
  
  # grab met
  met.subset <- met %>% 
    filter(date > forecast.start.day) %>% 
    filter(date <= forecast.end.day) 
  met.data <- met.subset %>% 
    select(all_of(met.driver))
  obs.temp <- met.subset %>% 
    select(min.temp.scale)
  gdd <- met.subset %>% 
    select(cum.gdd)
  
  n.days <- nrow(gdd) # number of days in forecast
  
  fore <- tick_forecast(params, ic, gdd, met.data, obs.temp, n.days) # run forecast
  predict <- fore[,check.day[t],] # pull out prediction for DA
  
  ### analysis step ###
  cat("\n-----------------------\n")
  index <- index + 1 # observation index counter
  
  loop.time <- Sys.time() - script.start
  cat("Elapsed time after", t, "forecasts\n")
  print(loop.time)
  cat(future.obs.index[t], "days forecasted \n")
  
  observation <- ticks.observed[,index] # ticks observed
  
  cat("\nLarva prediction summary:\n")
  print(summary(predict[1,]))
  cat("Larva observed", observation[1], "\n")
  cat("\nNymph prediction summary:\n")
  print(summary(predict[2,]))
  cat("Nymph observed", observation[2], "\n")
  cat("\nAdult prediction summary:\n")
  print(summary(predict[3,]))
  cat("Adult observed", observation[3], "\n\n")
  
  obs.prob <- obs_prob(ua, check.day[t], obs.temp[,1]) # observation probability for zero inflation
  obs.prob.mat <- matrix(NA, 3, Nmc) # storage 
  obs.prob.mat[1,] <- obs.prob$theta.larva[,check.day[t]]
  obs.prob.mat[2,] <- obs.prob$theta.nymph[,check.day[t]]
  obs.prob.mat[3,] <- obs.prob$theta.adult[,check.day[t]]
  
  # likelihood 
  cum.like <- rep(NA, Nmc) # store
  ens.like <- matrix(NA, Nmc, 3) # store
  for (i in 1:Nmc) {
    ens.like[i, ] <- dgpois(observation, predict[,i], 1-obs.prob.mat[,i]) # likelihood
    like.mu <- mean(ens.like[i, ]) # mean weight for each ensemble
    like[i, index] <- like.mu # store
    cum.like[i] <- prod(like[i, 1:index]) # convert to cumulative likelihood (weights)
  }
  
  # normalize and store weights to determine if we resample or not
  norm.weights <- cum.like / sum(cum.like)  # normalize weights
  like[, index] <- norm.weights             # store weights to carry over
  effect.size <- 1 / sum(norm.weights ^ 2)  # effective sample size
  cat("Effective Sample Size:", effect.size, "\n")
  
  # no resample
  if (effect.size > (Nmc / 2)) {
    cat("Resample: FALSE \n")
    resample[index] <- FALSE
    
    # weighted quantiles of ticks predicted
    wtd.quant <- matrix(NA, 3, 5)
    ic <- matrix(NA, Nmc, 3) # store
    for(i in 1:3){
      wtd.quant[i,] <- wtd.quantile(predict[i,], 
                                    cum.like / mean(cum.like),
                                    c(0.025,0.025,0.5,0.75,0.975))
      ic[,i] <- sample(wtd.quant[i,], Nmc, replace = TRUE)
    }
    ic[ic < 1] <- 0
    cat("Larva initial condition range:", range(ic[,1]), "\n")
    cat("Nymph initial condition range:", range(ic[,2]), "\n")
    cat("Adult initial condition range:", range(ic[,3]), "\n")
    
  } else { # resample
  
    cat("Resample: TRUE \n")
    resample[index] <- TRUE
    
    new.draw <- sample.int(Nmc, Nmc, replace = TRUE, prob = norm.weights) # resample with replacement
    params <- params[new.draw,] # new parameters
    params.hist[[index]] <- params # store parameters
    ic <- t(predict[,new.draw]) # new initial conditions
    ic[ic < 1] <- 0
    
    cat("Larva initial condition range:", range(ic[,1]), "\n")
    cat("Nymph initial condition range:", range(ic[,2]), "\n")
    cat("Adult initial condition range:", range(ic[,3]), "\n")
    
    like <- matrix(1, Nmc, length(days)) # reset weights
    weight.ncdf <- rep(1, Nmc) # weights for ncdf
    ua <- ua_parts(params, ic, c("parameter", "ic", "process")) # new uncertainty partitioning
  }
}

