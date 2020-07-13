library(LaplacesDemon)
library(plantecophys)
library(ecoforecastR)
library(tidyverse)
library(parallel)
library(abind)
library(mvtnorm)

script.start <- Sys.time()

source("Functions/site_data_met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/ua_parts.R")
source("Functions/obs_prob.R")
source("Functions/Future_Met.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/tick_forecast_parallel.R")
source("Functions/create_ncdf_tick.R")

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
Forecast_id <- "601144a5-ed57-4166-b943-38bfd66f9081"
ForecastProject_id <- 20200520 # Some ID that applies to a set of forecasts

## read array job number to paste into output file
grid.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# grid.id <- 1 # for test

# full site name, match to array job number
if(grid.id == 1) grid <- "Green Control" 
if(grid.id == 2) grid <- "Henry Control" 
if(grid.id == 3) grid <- "Tea Control" 

cat("\n--Running Tick Hindcast for", grid, "--\n")

grid.short <- gsub(" Control", "", grid) # site sub directory name
grid.no.space <- gsub(" ", "", grid) # site name without spaces

met.driver <- NULL

# load dir
in.dir <- "../FinalOut/A_Correct"
model.structure <- "RhoModels"
model.drivers <- "RhoAllStart" 
model.name <- paste0("Combined_thinMat_RhoAllStart_",
                     grid.no.space,
                     ".RData")

# load training fit
load(file.path(in.dir, 
               model.structure, 
               model.drivers, 
               grid.short, 
               model.name))

### parallel set-up
n.slots <- as.numeric(Sys.getenv("NSLOTS")) # read number of parallel slots (cores requested)
cat("Number of slots:", n.slots, "\n")
Nmc.per.node <- 5000
Nmc <- Nmc.per.node * n.slots

cat("Running", Nmc.per.node, "ensembles per node\n")
cat(Nmc, " total ensembles\n")

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

if(!dir.exists(out.dir)) dir.create(paste0(out.dir, "/"), recursive = TRUE)

# name for saving parameters
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

# t <- 1 # for testing
resample <- rep(NA, length(days))
weight.ncdf <- rep(1, Nmc) # weights for ncdf
norm.weights <- rep(log(1/Nmc), Nmc) # log scale
index <- 1
check.day <- rep(NA, length(days))
quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)

for (t in 1:(length(days)-1)) {
# for (t in 1:10) { # for testing
  
  ### forecast step ###
  forecast_issue_time <- as.Date(days[t])
  forecast.start.day <- as.character(days[t])  # date forecast issued
  forecast.end.day <- as.character(days[t+1]+16)  # next observation date
  check.day[t] <- as.integer(days[t+1] - days[t]) # day we evaluate forecast and run DA
  check.day[t] <- max(2, check.day[t]) # one-day forecasts evaluate 2 index as first is ic
  
  # always make at least 16-day forecast
  # if(check.day[t] < 16)  forecast.end.day <- as.character(days[t] + 16)
  
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
  
  # run forecast in parallel
  cl <- makeCluster(n.slots) # start cluster 
  datapieces <- clusterSplit(cl, 1:Nmc) # list to split data into equal chunks
  clusterExport(cl, c("params", "ic", "gdd", "met.data", "obs.temp", "n.days", "ua_parts", "obs_prob")) # export required objects
  fore <- parLapply(cl, datapieces, tick_forecast_parallel) # run forecast on cluster nodes
  stopCluster(cl) # stop cluster 
  
  # combine forecast output
  all.fore <- fore[[1]]
  for (slot in 2:n.slots) {
    all.fore <- abind(all.fore, fore[[slot]], along = 3)
  }
  predict <- all.fore[,check.day[t],] # pull out prediction for DA
  predict[predict == 0] <- 1E-10  # need to convert 0s for log.like calculation
  
  ### analysis step ###
  cat("\n-----------------------\n")
  index <- index + 1 # observation index counter
  
  loop.time <- Sys.time() - script.start
  cat("Elapsed time after", t, "forecasts\n")
  print(loop.time)
  cat(future.obs.index[t], "days forecasted \n")
  
  observation <- ticks.observed[,index] # ticks observed
  
  cat("\nLarva prediction summary:\n")
  print(round(quantile(predict[1,], quants)))
  cat("Larva observed", observation[1], "\n")
  cat("\nNymph prediction summary:\n")
  print(round(quantile(predict[2,], quants)))
  cat("Nymph observed", observation[2], "\n")
  cat("\nAdult prediction summary:\n")
  print(round(quantile(predict[3,], quants)))
  cat("Adult observed", observation[3], "\n\n")
  
  obs.prob <- obs_prob(ua, 1, obs.temp[check.day[t],1]) # observation probability for zero inflation
  obs.prob.mat <- matrix(NA, 3, Nmc) # storage 
  obs.prob.mat[1,] <- obs.prob$theta.larva[,1]
  obs.prob.mat[2,] <- obs.prob$theta.nymph[,1]
  obs.prob.mat[3,] <- obs.prob$theta.adult[,1]
  
  ## doing obs prob twice???
  
  # likelihood 
  cum.like <- rep(NA, Nmc) # store
  ens.like <- matrix(NA, Nmc, 3) # store
  for (i in 1:Nmc) {
    omega <- 1 - min(obs.prob.mat[,i], 1)
    ens.like[i, ] <- dgpois(observation, predict[,i], omega, log = TRUE) # likelihood
  }
  
  log.like.mu <- apply(ens.like, 1, mean) # mean weight for each ensemble
  w <- exp(norm.weights)*exp(log.like.mu) # previous normalized weights * current weights
  norm.weights <- w/sum(w) # normalize 

  print(sum(norm.weights))

  # normalize and store weights to determine if we resample or not
  # norm.weights <- cum.like / sum(cum.like)  # normalize weights
  # like[, index] <- norm.weights             # store weights to carry over
  effect.size <- 1 / sum(norm.weights ^ 2)  # effective sample size
  cat("Effective Sample Size:", effect.size, "\n")
  weight.ncdf <- norm.weights
  # logL = rep(NA, Nmc)
  # for (i in 1:Nmc) {
  #   logL[i] = mean(dgpois(observation, predict[,i], 1-obs.prob.mat[,i], log = T))
  # }
  # norm.weights = exp(logL) / sum(exp(logL))  ## weight each run
  # effect.size <- 1 / sum(norm.weights ^ 2)  # effective sample size
  # cat("Effective Sample Size:", effect.size, "\n")

  # no resample
  if (effect.size > (Nmc / 2)) {
    cat("Resample: FALSE \n")
    resample[index] <- FALSE
    
    # weighted quantiles of ticks predicted
    wtd.quant <- matrix(NA, 3, 5)
    ic <- matrix(NA, Nmc, 3) # store
    for(i in 1:3){
      wtd.quant[i,] <- wtd.quantile(predict[i,], 
                                    norm.weights / mean(norm.weights),
                                    c(0.025,0.025,0.5,0.75,0.975))
      ic[,i] <- sample(wtd.quant[i,], Nmc, replace = TRUE)
    }
    
    ic[ic < 1] <- 0
    cat("Larva initial condition summary:", 
        round(quantile(ic[,1], quants)), "\n")
    cat("Nymph initial condition range:", 
        round(quantile(ic[,2], quants)), "\n")
    cat("Adult initial condition range:",
        round(quantile(ic[,3], quants)), "\n")
    
  } else { # resample
  
    cat("Resample: TRUE \n")
    resample[index] <- TRUE
    
    new.draw <- sample.int(Nmc, Nmc, replace = TRUE, prob = norm.weights) # resample with replacement
    print(head(new.draw))
    params.new <- params[new.draw,-c(1:9)] # new parameters
    params.new.sigma <- params[new.draw,1:9] # grab prec matrix elements
    new.states <- t(predict[,new.draw]) # new initial conditions
    
    X <- cbind(new.states, params.new) # combine states and parameters
    Xbar <- apply(X, 2, mean)      # mean after resampling
    SIGMA <- cov(t(t(X) / Xbar))   # covariance matrix, normalized to the mean
    
    ## Kernel smoothing (to avoid ensemble members with identical parameters)
    h <- 0.9
    e <- rmvnorm(n = Nmc, mean = rep(0, ncol(X)), sigma = SIGMA)
    Xstar <- t(Xbar + h*(t(X)-Xbar) + t(e)*Xbar*(1-h))
    
    
    # ic after smoothing
    ic <- matrix(NA, Nmc, 3)
    for (i in 1:3) ic[, i] <- round(pmax(Xstar[, i], 0)) ## ticks are discrete and positive
    
    # ic <- new.states
    
    cat("Larva initial condition summary:", 
        round(quantile(ic[,1], quants)), "\n")
    cat("Nymph initial condition range:", 
        round(quantile(ic[,2], quants)), "\n")
    cat("Adult initial condition range:",
        round(quantile(ic[,3], quants)), "\n")
    
    # parameters after smoothing
    params <- cbind(params.new.sigma, Xstar[,-c(1:3)])
    params.hist[[index]] <- params # store parameters
    
    norm.weights <- rep(log(1/Nmc), Nmc) # reset weights
    ua <- ua_parts(params, ic, c("parameter", "ic", "process")) # new uncertainty partitioning
  }
  
  ncfname <- paste(t, "Tick_hindcast_pf.nc", sep = "_") # name file
  data_assimilation <- rep(0, n.days)
  data_assimilation[check.day[t]] <- 1
  
  # write netCDF
  create_ncdf_tick(
    file.path(out.dir, ncfname),
    all.fore,
    forecast_issue_time,
    n.days+1, # add one for IC storage
    Nmc,
    data_assimilation,
    weight.ncdf,
    ForecastProject_id,
    Forecast_id,
    forecast_issue_time)
}

# par(mfrow=c(2,3))
# hist(params.hist[[1]][,10])
# for(i in 2:6){
#   hist(params.hist[[i]][,10])
# }

save(params.hist, resample, check.day,
     file = file.path(out.dir, param.name))

cat("Parameters saved\n")
cat("--- END ---\n")