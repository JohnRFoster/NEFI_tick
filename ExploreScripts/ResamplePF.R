library(ecoforecastR)

## Load Model Function
source("Models/ensemble_gdd_k_window.R")
# source("Models/ensemble_temp_window_3alpha.R")
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

## out names
out.folder <- "../FinalOut/DA_Runs/K_estimate_4alpha/"

sites <- c("Green Control", "Henry Control", "Tea Control")
xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
# xx=1
site <- sites[xx]
cat("\n Resample PF for", site, "\n")
if(site == "Green Control"){
  s <- 1
} else if(site == "Henry Control"){
  s <- 2
} else {
  s <- 3
}
site.out <- gsub(" ", "", site)
# number of ensembles to run
Nmc <- 400 # 350 * 16 = 560 total ens members for each day

# run.type <- "Old"
run.type <- "New"

out.name <- paste(site.out, model.type, run.type, Nmc, sep = "_")

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
last.day <- as.Date(last.day[s,2])
future.days <- seq.Date(last.day + 1, last.day + (365*10), by = 1)

## index of when observations happen (days from last date)
obs.index <- vector()
vec <- which(ticks.future$date %in% future.days)
for(i in seq_along(which(ticks.future$date %in% future.days))){
  xx <- which(ticks.future$date[i] == future.days)
  obs.index <- c(obs.index, xx)
}

day <- 0
params.archive <- list()
resamp.archive <- list()
out.ci <- out.mean <- list()
resample.hist <- rep(NA, length(obs.index))
N <- Nmc/2

for(t in seq_along(future.days)){
# for(t in 1:50){
# for(t in 317:350){
  ## month counter
  if(t %% 31 == 0){print(t)}
  
  ## run ensemble
  gdd <- met$cum.gdd[t:(N_days + t - 1)]
  temp <- met$temp[t:(N_days + t - 1)]
  out  <- ensemble_gdd_k_window(site, params, ic, gdd, N_days, Nmc)
  # out  <- ensemble_temp_window_3alpha(site, params, ic, temp, gdd, N_days, Nmc)
  resamp.archive[[t]] <- out

  ## initial condition for next forecast
  if(t == 1){
    ic <- out[,1,]
    out.ci[[t]] <- t(apply(ic, 1, quantile, c(0.025, 0.5, 0.975)))
    out.mean[[t]] <- t(apply(ic, 1, mean))
    next
  }
  if(t <= 16){
    build <- t
    day.ens <- array(NA, dim = c(3, t, Nmc))
  } else {
    build <- N_days
    day.ens <- array(NA, dim = c(3, N_days, Nmc))
  }
  seq <- t+1-N_days:1
  # print(seq)
  seq <- seq[seq > 0]
  
  # loop over days that have forecasts
  for(d in seq){
    fore <- resamp.archive[[d]][,build,]
    for(ens in 1:Nmc){
      day.ens[,build,ens] <- fore[,ens] 
    }
    build <- build - 1
  }
  ic[1,] <- apply(day.ens[1,,], 2, median)
  ic[2,] <- apply(day.ens[2,,], 2, median)
  ic[3,] <- apply(day.ens[3,,], 2, median)
  out.ci[[t]] <- t(apply(day.ens, 1, quantile, c(0.025, 0.5, 0.975)))
  out.mean[[t]] <- t(apply(day.ens, 1, mean))

  ## analysis step - if observation is present
  if(t %in% obs.index){
    day <- day + 1
    cat("\nObservation", day, "\n")
    cat(day / length(obs.index) * 100, "% Done\n")
    
    ### old way
    # if(day == 1){
    #   wt <- rep(1, Nmc)
    # }
    # build <- N_days
    # day.ens <- array(NA, dim = c(3, N_days, Nmc))
    # 
    # # loop over days that have forecasts
    # for(d in seq_along(seq)){
    #   fore <- resamp.archive[[seq[d]]][,build,]
    #   for(ens in 1:Nmc){
    #     day.ens[,build,ens] <- fore[,ens]
    #   }
    #   build <- build - 1
    # }
    # 
    # like <- apply(day.ens, 3, median)
    # 
    # wt.new <- dpois(ticks.future$obs[, day], like)
    # wt <- wt*wt.new
    # effect.size <- 1 / sum(wt^2)
    # 
    # cat("Effective Size:", effect.size, "\n")
    # 
    # # ## resample ensembles
    # index <- sample.int(Nmc, Nmc, replace = TRUE, prob = wt)
    # 
    # obs.diff <- as.integer(obs.index[day] - obs.index[day-1])
    # if(day > 1){
    #   if(obs.diff > N_days){
    #     wt <- rep(1, Nmc)
    #   } else {
    #     wt.test <- apply(day.ens[,obs.diff:N_days,], 3, median)
    #   }
    # }
    # 
    
    ## think about cumulative weights through time
    ## resample when effective sample size drops below some threshold
    ## effective sample size = 1 / sum(wt^2)
    
    ## new way
    if(day == 1){
      wt <- rep(1, Nmc*N_days)
    }
    build <- N_days
    fore <- resamp.archive[[t-N_days]][,build,]
    for(d in t+1-N_days:2){
      build <- build-1
      fore <- cbind(fore, resamp.archive[[d]][,build,])
    }
    like <- apply(fore, 2, median)
    wt.new <- dpois(ticks.future$obs[, day], like)
    wt <- wt*wt.new

    ens.vec <- rep(1:Nmc, N_days)
    index <- base::sample(ens.vec, Nmc*N_days, replace = TRUE, prob = wt)
    rel.wt <- table(index) / (Nmc*N_days)
    rel.index <- rep(0, Nmc)
    rel.index[as.numeric(names(rel.wt))] <- as.numeric(rel.wt)
    index <- sample.int(Nmc, Nmc, replace = TRUE, prob = rel.index)
    # effect.size <- 1 / sum(wt^2)
    effect.size <- 1 / sum(rel.index^2)
    cat("Effective Size:", effect.size, "\n")
    
    if(day > 1){
      obs.diff <- obs.index[day] - obs.index[day-1]
      if(obs.diff > N_days){
        wt <- rep(1, Nmc) # reset wt
      } else {
        wt.keep <- obs.diff*Nmc
        wt[1:(wt.keep-1)] <- 1 # change what we discard to 1
      }
    }
  
   if(effect.size > N){ # dont resample
      resample.hist[day] <- 0
      params.archive[[day]] <- params
      ic[1,] <- apply(day.ens[1,,], 2, median)
      ic[2,] <- apply(day.ens[2,,], 2, median)
      ic[3,] <- apply(day.ens[3,,], 2, median)
      cat("Resample: FALSE \n")
    } else { # resample
      resample.hist[day] <- 1
      params <- params[index,]
      params.archive[[day]] <- params
      day.ens <- day.ens[,,index]
      ic[1,] <- apply(day.ens[1,,], 2, median)
      ic[2,] <- apply(day.ens[2,,], 2, median)
      ic[3,] <- apply(day.ens[3,,], 2, median)
      cat("Resample: True \n")
    }
  }
}

time.out <- as.Date(Sys.time(), "%Y%m%d")

obs <- ticks.future$obs
out.save <- paste0("Production", "CI.Mean.Median_", time.out, ".RData")
out.save <- paste(out.name, "N", N, out.save, sep = "_")
save(out.ci, out.mean, obs.index, obs, file = file.path(out.folder,out.save))

out.save <- paste0("Production", "resamp.archive_", time.out, ".RData")
out.save <- paste(out.name, "N", N, out.save, sep = "_")
save(resamp.archive, file = file.path(out.folder,out.save))

out.save <- paste0("Production", "params.archive_", time.out, ".RData")
out.save <- paste(out.name, "N", N, out.save, sep = "_")
save(params.archive, resample.hist, params.orig, file = file.path(out.folder,out.save))

