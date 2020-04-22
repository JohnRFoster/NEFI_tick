library(ecoforecastR)
library(boot)
library(lubridate)

source("Functions/mice_06_08.R")
source("Functions/Future_Met.R")
source("Functions/create_ncdf_mouse.R")
source("Functions/loglik_Bernoulli.R")
source("Functions/mouse_forecast.R")

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
Forecast_id <- "05cf9022-6ea7-415f-88b7-c4f546a596a0"
ForecastProject_id <-
  20200414 # Some ID that applies to a set of forecasts

dir <- "../FinalOut/DA_Runs/Mouse_Hindcast/Green"


new_mouse_observation <- function(capt.history, mice, days, index) {
  current.date <- days[index]
  
  if (current.date %in% mice$Full.Date.1) {
    current.matrix <- mice %>%
      filter(Full.Date.1 == current.date) %>%
      select(Tag.., Day.1) %>%
      filter(Day.1 == 1)
  } else if (current.date %in% mice$Full.Date.2) {
    current.matrix <- mice %>%
      filter(Full.Date.2 == current.date) %>%
      select(Tag.., Day.2)
  }
  
  # full capture history
  capt.history <-
    full_join(capt.history, current.matrix, by = "Tag..")
  
  # convert NAs to 0
  last <- ncol(capt.history)
  capt.history[which(is.na(capt.history[, last])), last] <- 0
  
  # capt.history <- capt.history[,c(1,ncol(capt.history))]
  
  return(capt.history)
}

load("../FinalOut/GreenControlMR/GreenControlMR_1_5.RData")
params <- as.matrix(jags.out)
Nmc <- 500
draw <- sample.int(nrow(params), Nmc)
params <- params[draw, ]

# store historical parameters
params.hist <- list()
params.hist[[1]] <- params


mice.future <- suppressWarnings(mice_06_18("Green Control"))
mice <- mice.future$table
ch <- mice.future$full.matrix
mice.observed <- colSums(ch)

day.1 <-
  as.character(unique(mice$Full.Date.1))        # 1st capture date of sampling occasion
day.2 <-
  as.character(unique(mice$Full.Date.2))        # 2nd capture date of sampling occasion
days <-
  ymd(c(rbind(day.1, day.2)))                    # vector of unique trapping days (for colnames)
future.obs.index <- cumsum(as.numeric(diff(days)))
future.obs.index <- c(1, future.obs.index + 1)

met <- future_met("Green Control")
met$date <- as.character(met$date)
start.date <- as.character(days[1])
start.index <- match(start.date, met$date)

met <- met[start.index:nrow(met), ]







## first initial conditions is first observation of 2006
current.date <- days[1]
current.matrix <- mice %>%
  filter(Full.Date.1 == current.date) %>%
  select(Tag.., Day.1) %>%
  filter(Day.1 == 1)

capt.history <- current.matrix

# first initialize every ensemble to the capture history
ic <- matrix(0, nrow(current.matrix) * 4, Nmc)
for (m in 1:Nmc) {
  ic.index <- sample(nrow(capt.history), sum(capt.history[, 2]))
  ic[ic.index, m] <- 1
  
  # total number of latent mice
  tot.mice <-
    round(sum(capt.history[, 2]) / params[m, "theta"]) - sum(capt.history[, 2])
  
  ic.index <- sample(nrow(capt.history):nrow(ic), tot.mice)
  ic[ic.index, m] <- 1
}


t <- 2 # for testing
n.days <- 16
index <- 1
weight <- rep(1, Nmc)
log.like <- matrix(0, Nmc, length(future.obs.index))
data_assimilation <- 0
store.forecast.Nmc <- ens.store <- list()
lambda <- matrix(1, Nmc, future.obs.index[length(future.obs.index)])
resample <- rep(NA, length(future.obs.index))
forecast.start.index <- 2
days.since.obs <- 0


for (t in 2:(365 * 2)) {
  # for(t in 2:23){
  
  ### forecast step ###
  # grab met
  met.seq <- forecast.start.index:(forecast.start.index + n.days - 1)
  precip <- met$precip[met.seq]
  temp <- met$temp.scale[met.seq]
  
  for (d in seq_along(temp)) {
    lambda[, met.seq[d]] <- inv.logit(params[, "lambda.mean"] +
                                        params[, "beta[1]"] * precip[d] +
                                        params[, "beta[2]"] * temp[d])
  }
  
  fore <- mouse_forecast(n.days, params, ic, precip, temp)
  
  # print(range(colSums(fore[2,,,n.days])))
  # first date of forecast
  forecast_issue_time <- as.Date(days[t - 1])
  
  # name file
  ncfname <- paste(t, "Mouse_hindcast_pf.nc", sep = "_")
  
  # write netCDF
  # source("Functions/create_ncdf_mouse.R")
  # create_ncdf_mouse(file.path(dir,ncfname),
  #                   fore,
  #                   forecast_issue_time,
  #                   n.days,
  #                   Nmc,
  #                   data_assimilation,
  #                   ForecastProject_id,
  #                   Forecast_id,
  #                   forecast_issue_time)
  
  # reset data assimilation to 0
  data_assimilation <- 0
  
  ### analysis step - if we make new observation ###
  if (t %in% future.obs.index) {
    cat("\n-----------------------\n")
    cat(t, "days forecasted \n")
    
    # assimilating data, set to 1
    # ??? might want to set to particle weights ???
    data_assimilation <- 1
    
    # observation index counter
    index <- index + 1
    
    # forecast start day counter
    forecast.start.index <- t + 1
    
    # observation
    capt.history <-
      new_mouse_observation(capt.history, mice, days, index)
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
      cat("Removing ",
          length(dead),
          "mice out of",
          nrow(capt.history.new),
          "in matrix\n")
      capt.history <- capt.history[-dead, ]
    }
    
    
    # most recent capture day only
    obs <- capt.history[, ncol(capt.history)]
    
    # total number observed
    obs.total <- sum(obs, na.rm = TRUE)
    
    # prediction
    if (days.since.obs > 2) {
      pred.full <- fore[1, , , days.since.obs]
    } else {
      pred.full <- fore[1, , , 2]
    }
    
    if (length(obs) < nrow(pred.full)) {
      pred <- pred.full[1:length(obs), ]
    }
    
    # total predicted
    pred.total <- colSums(fore[2, , , dim(fore)[4]])
    pred.quantile <- quantile(pred.total, c(0.025, 0.5, 0.975))
    
    ens.store[[index]] <- pred.total
    
    cat("\npredicted",
        pred.quantile,
        "\nobserved",
        obs.total,
        "\n")
    
    # likelihood (assign particles a weight)
    # like <- matrix(NA, length(obs), Nmc)
    # for(i in seq_along(obs)){ # individuals
    #   for(m in 1:Nmc){ # ensembles
    #     like[i,m] <- dbinom(obs[i], 1, pred[i,m]*params[m,"theta"])
    #   }
    # }
    
    # logLike <- rep(NA, Nmc)
    cum.like <- matrix(NA, Nmc, index)
    for (i in 1:Nmc) {
      # log likelihood
      log.like[i, index] <- loglik_Bernoulli(params[i, "theta"], pred[, i])
      
      # convert to cumulative likelihood
      cum.like[i, ] <- exp(cumsum(log.like[i, 1:index]))
    }
    
    weights <- cum.like[, index]
    
    # normalize weights
    norm.weights <- weights / sum(weights)
    print(round(norm.weights[1:20], 4))
    
    # effective sample size
    effect.size <- 1 / sum(norm.weights ^ 2)
    cat("Effective Sample Size:", effect.size, "\n")
    
    # reset
    n.days <- 16
    days.since.obs <- 0
    
    # no resample
    if (effect.size > (Nmc / 2)) {
      print(FALSE)
      resample[index] <- FALSE
      
      wts.ic <- matrix(NA, nrow(pred), 5)
      # ic <- matrix(0, length(obs)*4, Nmc)
      for (i in 1:nrow(pred)) {
        wts.ic <- wtd.quantile(pred.full[i, ],
                               weights / mean(weights),
                               c(0.025, 0.25, 0.5, 0.75, 0.975))
        
        top <- rep(c(wts.ic[1], wts.ic[2], wts.ic[3], wts.ic[4], wts.ic[5]),
              each = Nmc / 5)
        
        ic[i, ] <- top
      }
      cat("Initial Condition range:", range(colSums(ic)))
      
      # cat("Weighted quantiles:", wts.ic, "\n")
      
      # resample
    } else {
      print(TRUE)
      resample[index] <- TRUE
      
      # resample with replacement
      new.draw <- sample.int(Nmc, Nmc, replace = TRUE, prob = weight)
      
      # new parameters
      params <- params[new.draw, ]
      
      # store parameters
      params.hist[[index]] <- params
      
      # new initial conditions
      ic <- pred.full[, new.draw]
      
      cat("Initial Condition range:", range(colSums(ic)))
      
      # reset weights
      for (i in 1:Nmc) {
        log.like[i, 1:index] <- 0
      }
    }
  } else {
    # initial condition when no new observation
    # is the one-day ahead forecast
    # ic <- fore[1,,,2]
    
    
    # add a day to forecast
    n.days <- n.days + 1
    days.since.obs <- days.since.obs + 1
  }
}

par(mfrow = c(3, 3))
hist(params.hist[[1]][, "lambda.mean"])
for (i in 2:9) {
  hist(params.hist[[i]][, "lambda.mean"], main = i)
}

Nmc.normalized <- store.forecast.Nmc
quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)

diff.2 <- which(diff(future.obs.index) != 1)
one.day.index <- future.obs.index[diff.2]
two.day.index <- one.day.index + 1
other.index <- one.day.index - 1

grab.index <- sort(c(one.day.index, two.day.index, other.index))
grab.index <- grab.index[grab.index < 210]

plot.latent <- Nmc.normalized[[2]][1, , , 2]
plot.obs <- Nmc.normalized[[2]][2, , , 2]

day.sum <- quantile(colSums(plot.latent), quants)

for (d in 3:length(grab.index)) {
  grab <- grab.index[d]
  
  # only want second day of forecast
  if (grab %in% one.day.index) {
    plot.latent <- Nmc.normalized[[grab]][1, , , 2]
    plot.obs <- Nmc.normalized[[grab]][2, , , 2]
    
    day <- quantile(colSums(plot.latent), quants)
    day.sum <- rbind(day.sum, day)
    
  } else if (grab %in% (two.day.index)) {
    plot.latent <- Nmc.normalized[[grab]][1, , , 2]
    plot.obs <- Nmc.normalized[[grab]][2, , , 2]
    
    day <- quantile(colSums(plot.latent), quants)
    day.sum <- rbind(day.sum, day)
    
  } else {
    last.day <- dim(Nmc.normalized[[grab]])[4]
    plot.latent <- Nmc.normalized[[grab]][1, , , 2:last.day]
    plot.obs <- Nmc.normalized[[grab]][2, , , 2:last.day]
    
    
    for (t in 2:dim(plot.latent)[3]) {
      day <- quantile(colSums(plot.latent[, , t]), quants)
      day.sum <- rbind(day.sum, day)
    }
  }
}
str(day.sum)

mice.index <- future.obs.index[future.obs.index <= grab]
mice.2.plot <- mice.observed[1:length(mice.index)]

mice.vec <- rep(NA, nrow(day.sum))
mice.vec[mice.index] <- mice.2.plot

plot(1:nrow(day.sum),
     day.sum[, 5],
     pch = "",
     ylim = c(0, max(day.sum) + 5))
ciEnvelope(1:nrow(day.sum), day.sum[, 1], day.sum[, 5], col = "lightblue")
lines(1:nrow(day.sum), day.sum[, 3])
points(1:nrow(day.sum), mice.vec)


for (i in 1:length(params.hist)) {
  d <- quantile(params.hist[[i]][, "lambda.mean"])
  print(d)
}


boxplot(ens.store)

end <- length(ens.store)

points(1:end, mice.observed[1:end], pch = 16, col = "blue")

ens <- ens.store[[2]]
ens <- ens[ens >= quantile(ens, 0.025) & ens <= quantile(ens, 0.975)]

plot(2:end, pch = "", ylim = c(0, max(mice.observed) + 5))
points(rep(2, length(ens)), ens, pch = 18)

for (i in 3:length(ens.store)) {
  ens <- ens.store[[i]]
  ens <- ens[ens >= quantile(ens, 0.025) & ens <= quantile(ens, 0.975)]
  points(rep(i, length(ens)), ens, pch = 18)
}
points(2:end, mice.observed[2:36], pch = 16, col = "red")
