script.start <- Sys.time()

library(ecoforecastR)
library(boot)
library(lubridate)
library(parallel)
library(abind)

source("Functions/mice_06_08.R")
source("Functions/Future_Met.R")
source("Functions/create_ncdf_mouse.R")
source("Functions/mouse_forecast_parallel.R")
source("Functions/new_mouse_observation.R")

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
Forecast_id <- "05cf9022-6ea7-415f-88b7-c4f546a596a0"
ForecastProject_id <- 20200414 # Some ID that applies to a set of forecasts

dir <- "../FinalOut/DA_Runs/Mouse_Hindcast/Green"

file.name <- "Green_parameters.RData"
file <- file.path(dir, file.name)

print(file)

### parallel set-up
# read number of parallel slots (cores requested)
n.slots <- as.numeric(Sys.getenv("NSLOTS"))
cat("Number of slots:", n.slots, "\n")

# make clusters
# cl <- parallel::makeCluster(n.slots-1)
# showConnections()

load("../FinalOut/GreenControlMR/GreenControlMR_1_5.RData")
params <- as.matrix(jags.out)
Nmc.per.node <- 100
Nmc <- Nmc.per.node * n.slots
draw <- sample.int(nrow(params), Nmc)
params <- params[draw, ]

cat("Running", Nmc.per.node, "ensembles per node\n")
cat(Nmc, " total ensembles\n")

# store historical parameters
params.hist <- list()
params.hist[[1]] <- params


mice.future <- suppressWarnings(mice_06_18("Green Control"))
mice <- mice.future$table
ch <- mice.future$full.matrix
mice.observed <- colSums(ch)

day.1 <- as.character(unique(mice$Full.Date.1)) # 1st capture date of sampling occasion
day.2 <- as.character(unique(mice$Full.Date.2)) # 2nd capture date of sampling occasion
days <-  ymd(c(rbind(day.1, day.2))) # vector of unique trapping days (for colnames)
future.obs.index <- cumsum(as.numeric(diff(days)))
future.obs.index <- c(1, future.obs.index + 1)

met <- future_met("Green Control")
met$date <- as.character(met$date)
start.date <- as.character(days[1])
start.index <- match(start.date, met$date)

met <- met[start.index:nrow(met), ]

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
n.days <- 16
index <- 1
weights <- rep(1, Nmc)
like <- matrix(1, Nmc, length(future.obs.index))
data_assimilation <- 0
store.forecast.Nmc <- ens.store <- all.forecasts <- list()
lambda <- matrix(1, Nmc, future.obs.index[length(future.obs.index)])
resample <- rep(NA, length(future.obs.index))
forecast.start.index <- 2
days.since.obs <- 0
weight.ncdf <- rep(1, Nmc) # weights for ncdf



for (t in 2:(365*10)) {
  # for(t in 2:126){
  
  ### forecast step ###
  
  # grab met
  met.seq <- forecast.start.index:(forecast.start.index + n.days - 1)
  precip <- met$precip[met.seq]
  temp <- met$temp.scale[met.seq]
  
  # calculate and store daily survival
  for (d in seq_along(temp)) {
    lambda[, met.seq[d]] <- inv.logit(params[, "lambda.mean"] +
                                        params[, "beta[1]"] * precip[d] +
                                        params[, "beta[2]"] * temp[d])
  }
  # start cluster (move to beginning of script???)
  cl <- makeCluster(n.slots)
  
  # list to split data into equal chunks
  datapieces <- clusterSplit(cl, 1:Nmc)
  
  # export required objects
  clusterExport(cl, c("n.days", "params", "ic", "precip", "temp"))
  
  # run forecast on cluster nodes
  fore <- parLapply(cl, datapieces, mouse_forecast_parallel)
  
  # stop cluster (move to end of script???)
  stopCluster(cl)
  
  # combine forecast output
  all.fore <- fore[[1]]
  for (slot in 2:n.slots){
    all.fore <- abind(all.fore, fore[[slot]], along = 3)
  }

  # first date of forecast
  forecast_issue_time <- as.Date(days[t - 1])
  
  # reset data assimilation to 0
  data_assimilation <- 0
  
  ### analysis step - if we make new observation ###
  if (t %in% future.obs.index) {
    
    all.forecasts[[index]] <- fore
    
    cat("\n-----------------------\n")
    loop.time <- Sys.time() - script.start
    cat("Elapsed time:\n")
    print(loop.time)
    cat(t, "days forecasted \n")
    
    # assimilating data, set to 1
    # ??? might want to set to particle weights ???
    data_assimilation <- 1
    
    # observation index counter
    index <- index + 1
    
    # forecast start day counter
    forecast.start.index <- t + 1
    
    # observation
    capt.history <- new_mouse_observation(capt.history, mice, days, index)
    capt.history.new <- capt.history
    
    # find last day each mouse was captured
    dead <- 1
    for (r in 1:nrow(capt.history.new)) {
      check <- suppressWarnings(max(which(capt.history.new[r, ] == 1)))
      capture.day <- future.obs.index[check]
      if ((t - capture.day) > 50 & !is.infinite(check)) {
        # cat("Checking survival for mouse", capt.history[r,"Tag.."], "\n")
        surv.prob <- rep(NA, Nmc)
        for (m in 1:Nmc) {
          surv.prob[m] <- prod(lambda[m, check:t])
        }
        if (median(surv.prob) < 0.05) {
          dead <- c(dead, r)
        }
      }
    }
    if (length(dead) > 1) {
      dead <- dead[-1]
      cat("Removing ", length(dead), "mice out of", nrow(capt.history.new), "in matrix\n")
      capt.history <- capt.history[-dead, ]
    }
    
    
    # most recent capture day only
    obs <- capt.history[, ncol(capt.history)]
    
    # total number observed
    obs.total <- sum(obs, na.rm = TRUE)
    
    # prediction and store
    if (days.since.obs > 2) {
      store.forecast.Nmc[[index]] <- all.fore
      pred.full <- all.fore[1, , , days.since.obs]
      pred.obs <- all.fore[2, , , days.since.obs]
    } else {
      store.forecast.Nmc[[index]] <- all.fore[,,,2]
      pred.full <- all.fore[1, , , 2]
      pred.obs <- all.fore[2, , , 2]
    }
    
    if (length(obs) < nrow(pred.full)) {
      pred <- pred.obs[1:length(obs), ]
      pred.latent <- pred.full[1:length(obs), ]
    }
    
    # total predicted
    pred.total <- colSums(pred.obs)
    pred.quantile <- quantile(pred.total, c(0.025, 0.5, 0.975))
    
    ens.store[[index]] <- pred.total
    
    cat("\npredicted summary:\n")
    print(summary(pred.total))
    cat("observed", obs.total, "\n")
    
    # likelihood (assign particles a weight)
    cum.like <- rep(NA, Nmc)
    ens.like <- matrix(NA, Nmc, length(obs))
    for (i in 1:Nmc) {
      
      # likelihood
      ens.like[i,] <- dbinom(obs, 1, params[i, "theta"] * pred.latent[, i])
      
      # mean weight for each ensemble
      like.mu <- mean(ens.like[i,])
      
      # store
      like[i, index] <- like.mu
      
      # convert to cumulative likelihood (weights)
      cum.like[i] <- prod(like[i, 1:index])
    }
    
    # normalize weights
    norm.weights <- cum.like / sum(cum.like)
    
    # store weights to carry over
    like[,index] <- norm.weights
    
    # effective sample size
    effect.size <- 1 / sum(norm.weights ^ 2)
    cat("Effective Sample Size:", effect.size, "\n")
    
    # reset
    n.days <- 16
    days.since.obs <- 0
    
    # no resample
    if (effect.size > (Nmc / 2)) {
      cat("Resample: FALSE \n")
      resample[index] <- FALSE
      
      # mice we know are alive
      mice.1 <- which(obs == 1)
      
      # number of mice to inflate
      # inflate <- round(pred.total / params[, "theta"])
      
      # weighted quantiles of total latent mice predicted
      wts.ic <- wtd.quantile(pred.total,
                             weights / mean(weights),
                             c(0.025, 0.25, 0.5, 0.75, 0.975))
      
      build.mice <- rep(wts.ic, each = Nmc / 5)
       
      ic <- matrix(0, length(obs)*4, Nmc)
      for (i in 1:Nmc) {
        ic[sample(nrow(ic), build.mice[i]), i] <- 1
        ic[mice.1, i] <- 1
      }
      cat("Initial Condition range:", range(colSums(ic)), "\n")
      cat("Dim IC:", dim(ic), "\n")
      
      # weights for ncdf
      weight.ncdf <- norm.weights
      
    } else { # resample
      cat("Resample: TRUE \n")
      resample[index] <- TRUE
      
      # resample with replacement
      new.draw <- sample.int(Nmc, Nmc, replace = TRUE, prob = weights)
      
      # new parameters
      params <- params[new.draw, ]
      
      # store parameters
      params.hist[[index]] <- params
      
      # new predicted totals
      pred.total.draw <- pred.total[new.draw]
      
      # mice we know are alive
      mice.1 <- which(obs == 1)
      
      # new initial conditions
      ic <- matrix(0, length(obs)*4, Nmc)
      for(i in 1:Nmc){
        ic[1:length(obs),] <- obs
        
        # inflate to latent number
        inflate <- round(pred.total[new.draw[i]] / params[i, "theta"]) - pred.total[new.draw[i]]
        
        new.range <- (length(obs) + 1):nrow(ic) 
        ic[sample(new.range, inflate), i] <- 1 
        ic[mice.1, i] <- 1
      }
      
      cat("Initial Condition range:", range(colSums(ic)), "\n")
      cat("Dim IC:", dim(ic), "\n")
      
      # reset weights
      like <- matrix(1, Nmc, length(future.obs.index))
      
      # weights for ncdf
      weight.ncdf <- rep(1, Nmc)
    }
  } else {
    
    n.days <- n.days + 1 # add a day to forecast
    days.since.obs <- days.since.obs + 1 # add days since last observation
  }
  
  # name file
  ncfname <- paste(t, "Mouse_hindcast_pf.nc", sep = "_")
  
  # cat("Forecast dimensions", dim(all.fore), "\n")
  # cat("Days", n.days, "\n")
  # cat("Ens", Nmc, "\n")
  
  # write netCDF
  source("Functions/create_ncdf_mouse.R")
  create_ncdf_mouse(file.path(dir,ncfname),
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

save(params.hist, resample, future.obs.index,
     file = file)

cat("Cluster stopped\n")
cat("--- END ---\n")

# stop clusters
# parallel::stopCluster(cl)
# showConnections()


