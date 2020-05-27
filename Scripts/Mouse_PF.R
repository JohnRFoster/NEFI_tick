script.start <- Sys.time()

library(ecoforecastR)
library(boot)
library(lubridate)
library(parallel)
library(abind)
library(mvtnorm)

source("Functions/mice_06_08.R")
source("Functions/Future_Met.R")
source("Functions/create_ncdf_mouse.R")
source("Functions/mouse_forecast_parallel.R")
source("Functions/new_mouse_observation.R")

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
Forecast_id <- "05cf9022-6ea7-415f-88b7-c4f546a596a0"
ForecastProject_id <- 20200519 # Some ID that applies to a set of forecasts

grid <- "Green Control"
grid.short <- gsub(" Control", "", grid)

dir <- file.path("../FinalOut/DA_Runs/Mouse_Hindcast", grid.short)

file.name <- paste0(grid.short, "_parameters.RData")
file <- file.path(dir, file.name)

cat("Saving parameters to:\n", file, "\n")

### parallel set-up
n.slots <- as.numeric(Sys.getenv("NSLOTS")) # read number of parallel slots (cores requested)
cat("Number of slots:", n.slots, "\n")

# load parameters
if(grid == "Green Control") load("../FinalOut/GreenControlMR/GreenControlMR_1_5.RData")  
if(grid == "Henry Control") load("../FinalOut/HenryControlMR/HenryControlMR_1_6.RData")  
if(grid == "Tea Control") load("../FinalOut/TeaControlMR/TeaControlMR_1_6.RData")

params <- as.matrix(jags.out)
Nmc.per.node <- 125
Nmc <- Nmc.per.node * n.slots
draw <- sample.int(nrow(params), Nmc)
params <- params[draw, ]

cat("Running", Nmc.per.node, "ensembles per node\n")
cat(Nmc, " total ensembles\n")

# store historical parameters
params.hist <- list()
params.hist[[1]] <- params

mice.future <- suppressWarnings(mice_06_18(grid))
mice <- mice.future$table
ch <- mice.future$full.matrix
mice.observed <- colSums(ch)

day.1 <- as.character(unique(mice$Full.Date.1)) # 1st capture date of sampling occasion
day.2 <- as.character(unique(mice$Full.Date.2)) # 2nd capture date of sampling occasion
days <-  ymd(c(rbind(day.1, day.2))) # vector of unique trapping days (for colnames)
days <- days[days <= "2017-12-31"] # run hindcast through 2017 only, 2018 met data is NA
diff.days <- diff(days)
forecast.time.step <- pmax(16, diff(days)) # days to forecast

future.obs.index <- cumsum(as.numeric(diff(days)))
future.obs.index <- c(1, future.obs.index + 1)

# load and subset future met
met <- future_met(grid)
met$date <- as.character(met$date)
start.date <- as.character(days[1])
start.index <- match(start.date, met$date)
met <- met[start.index:nrow(met), ]

# calculate and store daily survival
lambda <- matrix(NA, Nmc, nrow(met))
for (d in 1:Nmc) {
  lambda[d,] <- inv.logit(params[d, "lambda.mean"] +
                          params[d, "beta[1]"] * met$precip +
                          params[d, "beta[2]"] * met$max.temp.scale)
}



## first initial conditions is first observation of 2006
current.date <- days[1]
capt.history <- mice %>%
  filter(Full.Date.1 == current.date) %>%
  select(Tag.., Day.1) %>%
  filter(Day.1 == 1)

# just want observation vector
obs <- capt.history[,2]

# first initialize every ensemble to the capture history
ic <- matrix(0, length(obs) * 4, Nmc)
for(i in 1:Nmc){
  ic[1:length(obs),] <- obs
  
  # inflate to latent number
  inflate <- round(sum(obs) / params[i, "theta"]) - sum(obs)
  
  new.range <- (length(obs) + 1):nrow(ic) 
  ic[sample(new.range, inflate), i] <- 1 
}

cat("Initial Condition range:", range(colSums(ic)), "\n")

t <- 2 # for testing
like <- matrix(1, Nmc, length(future.obs.index))
resample <- rep(NA, length(future.obs.index))
weight.ncdf <- rep(1, Nmc) # weights for ncdf
index <- 1
check.day <- rep(NA, length(days))

for (t in 1:(length(days)-1)) {

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
  precip <- met.subset %>% 
    select(precip)
  temp <- met.subset %>%
    select(max.temp.scale)
  
  n.days <- nrow(temp) # number of days in forecast
  
  # run forecast in parallel
  cl <- makeCluster(n.slots) # start cluster 
  datapieces <- clusterSplit(cl, 1:Nmc) # list to split data into equal chunks
  clusterExport(cl, c("n.days", "params", "ic", "precip", "temp")) # export required objects
  fore <- parLapply(cl, datapieces, mouse_forecast_parallel) # run forecast on cluster nodes
  stopCluster(cl) # stop cluster 
  
  # combine forecast output
  all.fore <- fore[[1]]
  for (slot in 2:n.slots) {
    all.fore <- abind(all.fore, fore[[slot]], along = 3)
  }

  ### analysis step ###
  cat("\n-----------------------\n")
  index <- index + 1 # observation index counter
  
  loop.time <- Sys.time() - script.start
  cat("Elapsed time after", t, "forecasts\n")
  print(loop.time)
  cat(future.obs.index[t], "days forecasted \n")
  
  # observation
  capt.history <- new_mouse_observation(capt.history, mice, days, index)
  capt.history.new <- capt.history
  
  # find last day each mouse was captured
  dead <- 1
  for (r in 1:nrow(capt.history.new)) {
    last.observed <- suppressWarnings(max(which(capt.history.new[r,-1] == 1)))
    capture.day <- future.obs.index[last.observed]
    if (!is.infinite(last.observed)) {
      surv.prob <- rep(NA, Nmc)
      for (m in 1:Nmc) {
        surv.prob[m] <- prod(lambda[m, capture.day:future.obs.index[t+1]])
      }
      if (median(surv.prob) < 0.05) {
        dead <- c(dead, r)
      }
    }
  }
  if (length(dead) > 1) {
    dead <- dead[-1]
    cat("Removing ", length(dead), "mice out of", nrow(capt.history.new), "in matrix\n")
    capt.history <- capt.history[-dead,]
  }
  
  obs <- capt.history[, ncol(capt.history)] # most recent capture day only
  obs.total <- sum(obs, na.rm = TRUE)  # total number observed
  
  # subset to correct day
  pred.full <- all.fore[1, , , check.day[t]] # predicted latent
  pred.obs <- all.fore[2, , , check.day[t]]  # predicted observed
  
  # remove annoted mice
  if (length(obs) < nrow(pred.full)) {
    pred <- pred.obs[1:length(obs),]
    pred.latent <- pred.full[1:length(obs),]
  }
  
  pred.total <- colSums(pred.obs) # total predicted
  pred.quantile <- quantile(pred.total, c(0.025, 0.5, 0.975)) # quantiles 
  
  cat("\npredicted summary:\n")
  print(summary(pred.total))
  cat("observed", obs.total, "\n")
  
  # likelihood 
  cum.like <- rep(NA, Nmc) # store
  ens.like <- matrix(NA, Nmc, length(obs)) # store
  for (i in 1:Nmc) {
    ens.like[i, ] <- dbinom(obs, 1, params[i, "theta"] * pred.latent[, i]) # likelihood
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
    
    mice.1 <- which(obs == 1) # mice we know are alive
    
    # weighted quantiles of total latent mice predicted
    wts.ic <- wtd.quantile(pred.total,
                           cum.like / mean(cum.like),
                           c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    build.mice <- rep(wts.ic, each = Nmc / 5) # helper for building ic matrix
    
    # build ic
    ic <- matrix(0, length(obs) * 4, Nmc)
    for (i in 1:Nmc) {
      ic[sample(nrow(ic), build.mice[i]), i] <- 1
      ic[mice.1, i] <- 1
    }
    cat("Initial Condition range:", range(colSums(ic)), "\n")
    cat("Dim IC:", dim(ic), "\n")
    
    weight.ncdf <- norm.weights # weights for ncdf
    
  } else { # resample
    
    cat("Resample: TRUE \n")
    resample[index] <- TRUE
    
    new.draw <- sample.int(Nmc, Nmc, replace = TRUE, prob = norm.weights) # resample with replacement
    params <- params[new.draw,] # new parameters
    pred.total.draw <- pred.total[new.draw] # new predicted totals
    
    ## Kernel smoothing (to avoid ensemble members with identical parameters)
    X <- cbind(pred.total.draw, params) # combine states and parameters
    Xbar <- apply(X, 2, mean)      # mean after resampling
    SIGMA <- cov(t(t(X) / Xbar))   # covariance matrix, normalized to the mean
    h <- 0.9
    e <- rmvnorm(n = Nmc, mean = rep(0, ncol(X)), sigma = SIGMA)
    Xstar = t(Xbar + h*(t(X)-Xbar) + t(e)*Xbar*(1-h))
    
    params <- Xstar[,-1]
    params.hist[[index]] <- params # store parameters
    
    mice.1 <- which(obs == 1) # mice we know are alive
    
    # new initial conditions
    ic <- matrix(0, length(obs) * 4, Nmc)
    for (i in 1:Nmc) {
      ic[1:length(obs), ] <- obs
      
      # inflate to latent number
      inflate <- round(pred.total[new.draw[i]] / params[i, "theta"]) - pred.total[new.draw[i]]
      
      new.range <- (length(obs) + 1):nrow(ic)
      ic[sample(new.range, inflate), i] <- 1
      ic[mice.1, i] <- 1
    }
    
    cat("Initial Condition range:", range(colSums(ic)), "\n")
    cat("Dim IC:", dim(ic), "\n")
    
    like <- matrix(1, Nmc, length(future.obs.index)) # reset weights
    weight.ncdf <- rep(1, Nmc) # weights for ncdf
    
    # re-calculate and store daily survival
    lambda <- matrix(NA, Nmc, nrow(met))
    for (d in 1:Nmc) {
      lambda[d,] <- inv.logit(params[d, "lambda.mean"] +
                                params[d, "beta[1]"] * met$precip +
                                params[d, "beta[2]"] * met$max.temp.scale)
    }
  }
  
  ncfname <- paste(t, "Mouse_hindcast_pf.nc", sep = "_") # name file
  data_assimilation <- rep(0, n.days)
  data_assimilation[check.day[t]] <- 1
  
  # write netCDF
  create_ncdf_mouse(
    file.path(dir, ncfname),
    all.fore,
    forecast_issue_time,
    dim(all.fore)[4],
    Nmc,
    data_assimilation,
    weight.ncdf,
    ForecastProject_id,
    Forecast_id,
    forecast_issue_time)
}

save(params.hist, resample, check.day,
     file = file)

cat("Parameters saved\n")
cat("--- END ---\n")

# stop clusters
# parallel::stopCluster(cl)
# showConnections()


