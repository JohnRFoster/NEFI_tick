library(ecoforecastR)
library(tidyverse)
library(plantecophys)

source("Functions/ua_parts.R")
source("Functions/obs_prob.R")
source("Functions/Future_Met.R")
source("Functions/get_ticks_2006_2018.R")
source("Models/jagsFilter_monthEffectLifeStage.R")
source("Functions/convergence_check.R")

script.start <- Sys.time()

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
# Forecast_id <- "601144a5-ed57-4166-b943-38bfd66f9081"
# ForecastProject_id <- 20200520 # Some ID that applies to a set of forecasts

## read array job number to paste into output file
grid.id <- 2 # for test
grid.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# full site name, match to array job number
grids <- c("Green Control", "Henry Control", "Tea Control")
grid <- grids[grid.id]

cat("\n-- Running Tick Hindcast for", grid, "--\n")

grid.short <- gsub(" Control", "", grid) # site sub directory name
grid.no.space <- gsub(" ", "", grid) # site name without spaces

met.driver <- NULL

# load dir
in.dir <- "../FinalOut/A_Correct"
model.structure <- "RhoModels"
model.drivers <- "RhoAllStart/MonthEffect" 
model.name <- paste0("Combined_thinMat_RhoStartMonthEffectLifeStage_",
                     grid.no.space,
                     ".RData")

# load training fit
load(file.path(in.dir, 
               model.structure, 
               model.drivers, 
               grid.short, 
               model.name))

# output directory
specific.name <- "LifeStageMonthEffect_UpdateAllParameters"
out.dir <- file.path("../FinalOut/DA_Runs/Tick_Hindcast/jagsFilter", 
                     model.structure,
                     model.drivers,
                     specific.name,
                     grid.short)

if(!dir.exists(out.dir)) dir.create(paste0(out.dir, "/"), recursive = TRUE)

# known <- "SIGMA" # currently only set up for SIGMA, but this setup should extend to other params
known <- NULL
if(!is.null(known)){ 
  known.vals <- list() # known values for jags data list
  if("SIGMA" %in% known){
    OMEGA <- params.mat[,grep(known, colnames(params.mat))]
    q.bar <- matrix(apply(OMEGA, 2, mean), 3, 3)  # Mean Omega, Precision 
    known.vals$SIGMA <- q.bar
  }
}


## Load Future ticks
ticks <- get_ticks_2006_2018(grid)
days <- ticks$date # date of observations
days <- days[days <= "2017-12-31"] # run hindcast through 2017 only, 2018 met data is NA
ticks.observed <- ticks$obs[,1:length(days)] # observed ticks

# load and subset future met
met <- future_met(grid)
met$date <- as.character(met$date)
start.date <- as.character(days[1])
start.index <- match(start.date, met$date)
met <- met[start.index:nrow(met), ]
met <- met %>% 
  select(c(date, cum.gdd, min.temp.scale)) # for katie csv
# met$vpd <- RHtoVPD(met$min.rh, met$min.temp)

# data and initial conditions for first hindcast
data <- update_data(params.mat, 4:12)
data$l.ic <- round(rpois(1, ticks.observed[1,1])) + 1
data$n.ic <- round(rpois(1, ticks.observed[2,1])) + 1
data$a.ic <- round(rpois(1, ticks.observed[3,1])) + 1

t <- 1
n.adapt <- 2000
n.iter <- 10000
n.chains <- 3
iter2save <- 10000
end <- (length(days)-1)
# end <- 3 # for testing

for (t in 1:end) {

  cat("==================================\n")
  cat(round(t/end*100, 2), "% done\n")
  
  ### forecast step ###
  forecast_issue_time <- as.Date(days[t])
  forecast.start.day <- as.character(days[t])  # date forecast issued
  forecast.end.day <- as.character(days[t+1])  # next observation date
  check.day <- as.integer(days[t+1] - days[t]) # day we evaluate forecast and run DA

  all.days <- seq.Date(ymd(forecast.start.day)+1, ymd(forecast.end.day), 1)
  data$month.index <- month(all.days)
  
  # grab met
  met.subset <- met %>% 
    filter(date > forecast.start.day) %>% 
    filter(date <= forecast.end.day) 
  obs.temp <- met.subset %>% 
    select(min.temp.scale)
  gdd <- met.subset %>% 
    select(cum.gdd)
  
  data$n.days <- nrow(gdd) # number of days in forecast
  data$gdd <- gdd
  data$met.obs <- obs.temp

  y <- matrix(NA, 3, check.day)
  y[,1] <- ticks.observed[,t]
  data$y <- y
  data$seq.days <- (data$n.days-1):1

  # run filter
  jags.filter <- run_jagsFilter(data, n.adapt, n.chains, known.vals)

  jags.out <- coda.samples(model = jags.filter$j.model,
                           variable.names = jags.filter$monitor,
                           n.iter = n.iter)

  cat("coda samples done, checking mcmc \n")
  
  ## convergence check
  out <- convergence_check(jags.out, jags.filter)
  
  preds <- as.matrix(out$predict)
  larva <- preds[,seq(1, by = 3, length.out = ncol(preds)/3)]
  nymph <- preds[,seq(2, by = 3, length.out = ncol(preds)/3)]
  adult <- preds[,seq(3, by = 3, length.out = ncol(preds)/3)]
  
  # quantiles for cat statements
  l.q <- round(quantile(larva[,ncol(larva)], c(0.025, 0.5, 0.975)))
  n.q <- round(quantile(nymph[,ncol(nymph)], c(0.025, 0.5, 0.975)))
  a.q <- round(quantile(adult[,ncol(adult)], c(0.025, 0.5, 0.975)))
  
  # cat statements for batch jobs
  cat("Larva oberved:", ticks.observed[1,t+1], "\n")
  cat("Larva forecast CI:", l.q, "\n\n")
  cat("Nymph oberved:", ticks.observed[2,t+1], "\n")
  cat("Nymph forecast CI:", n.q, "\n\n")
  cat("Adult oberved:", ticks.observed[3,t+1], "\n")
  cat("Adult forecast CI:", a.q, "\n\n")
  
  parameters <- as.matrix(out$params)
  data <- update_data(parameters, 1:12) # update with full posterior
  data$l.ic <- round(mean(larva[,ncol(larva)]))
  data$n.ic <- round(mean(nymph[,ncol(nymph)]))
  data$a.ic <- round(mean(adult[,ncol(adult)]))
  
  # thin for saving
  thin <- seq(1, nrow(parameters), length.out = iter2save)
  preds <- preds[thin,]
  parameters <- parameters[thin,]
  
  # save as .RData
  outname <- paste(t, specific.name, "Tick_hindcast_jagsFilter.RData", sep = "_") # name file
  save(preds, parameters, data,
       file = file.path(out.dir, outname))
    
}
