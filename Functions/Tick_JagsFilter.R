library(ecoforecastR)
library(tidyverse)
library(plantecophys)

source("Functions/ua_parts.R")
source("Functions/obs_prob.R")
source("Functions/Future_Met.R")
source("Functions/get_ticks_2006_2018.R")
source("Functions/run_enkf.R")
source("Functions/convergence_check.R")

script.start <- Sys.time()

# Forecast_id <- uuid::UUIDgenerate() # ID that applies to the specific forecast
# Forecast_id <- "601144a5-ed57-4166-b943-38bfd66f9081"
# ForecastProject_id <- 20200520 # Some ID that applies to a set of forecasts

## read array job number to paste into output file
# grid.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
grid.id <- 1 # for test

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

# Nmc <- nrow(params.mat)

update_data <- function(params.mat, predict.mat){
  
  # moment match beta distributions
  beta_moment_match <- function(param){
    ex <- mean(param)
    v <- var(param)
    alpha <- (((ex^2) * (1-ex)) / v) - ex
    beta <- (((ex*(1-ex)) / v) - 1) * (1-ex)
    moments <- list(alpha = alpha, beta = beta)
    return(moments)
  }
  
  data <- list() # data list
  prior.theta.adult <- beta_moment_match(params.mat[,"theta.adult"])
  prior.theta.nymph <- beta_moment_match(params.mat[,"theta.nymph"])
  data$beta.adult.alpha <- prior.theta.adult$alpha
  data$beta.adult.beta <- prior.theta.adult$beta
  data$beta.nymph.alpha <- prior.theta.nymph$alpha
  data$beta.nymph.beta <- prior.theta.nymph$beta
  data$beta.l.obs.mu <- mean(params.mat[,"beta.l.obs"])
  data$beta.l.obs.prec <- 1 / var(params.mat[,"beta.l.obs"])
  data$larva.mu <- mean(params.mat[,"phi.l.mu"])
  data$larva.prec <- 1 / var(params.mat[,"phi.l.mu"])
  data$nymph.mu <- mean(params.mat[,"phi.n.mu"])
  data$nymph.prec <- 1 / var(params.mat[,"phi.n.mu"])
  data$adult.mu <- mean(params.mat[,"phi.a.mu"])
  data$adult.prec <- 1 / var(params.mat[,"phi.a.mu"])
  data$l2n.mu <- mean(params.mat[,"grow.ln.mu"])
  data$l2n.prec <- 1 / var(params.mat[,"grow.ln.mu"])
  data$n2a.mu <- mean(params.mat[,"grow.na.mu"])
  data$n2a.prec <- 1 / var(params.mat[,"grow.na.mu"])
  data$a2l.mu <- mean(params.mat[,"repro.mu"])
  data$a2l.prec <- 1 / var(params.mat[,"repro.mu"])
  data$rho.l.mu <- mean(params.mat[,"rho.l"])
  data$rho.l.prec <- 1 / var(params.mat[,"rho.l"])
  data$rho.n.mu <- mean(params.mat[,"rho.n"])
  data$rho.n.prec <- 1 / var(params.mat[,"rho.n"])
  data$rho.a.mu <- mean(params.mat[,"rho.a"])
  data$rho.a.prec <- 1 / var(params.mat[,"rho.a"])  
  
  ## need to decompose SIGMA?
  
  # omega[t+1] <- q.bar[t] * beta[t+1]
  # beta[t+1] <- mean((omega[r,c]^2 + omega[r,r]*omega[c,c]) / var(omega[r,c]))
  
  # omega: process precision to the process covariance Q[t]
  # r and c: rows and columns of the process precision matrix
  # q.bar: posterior mean estimate of Q[t]
  
  wish_df <- function(Om, X, i, j, col) {
    (Om[i, j]^2 + Om[i, i] * Om[j, j]) / var(X[, col])
  }
  
  # process error draws, precision
  pattern <- "SIGMA"
  OMEGA <- params.mat[,grep(pattern, colnames(params.mat))]
  
  q.bar <- matrix(apply(OMEGA, 2, mean), 3, 3)  # Mean Omega, Precision
  
  col <- matrix(1:3 ^ 2, 3, 3)
  WV <- matrix(NA, 3, 3)
  for(r in 1:3){
    for(c in 1:3){
      WV[r, c] <- wish_df(q.bar, OMEGA, r, c, col[r, c])
    }
  }
  
  beta <- mean(WV)
  V <- solve(q.bar) * beta 
  
  # deal with posterior of sigma later
  data$R <- V
  data$k <- beta
  return(data)
}

# output directory
out.dir <- file.path("../FinalOut/DA_Runs/Tick_Hindcast/EnKF", 
                     model.structure,
                     model.drivers,
                     grid.short)

if(!dir.exists(out.dir)) dir.create(paste0(out.dir, "/"), recursive = TRUE)

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


data <- update_data(params.mat, predict.mat)
data$l.ic <- round(rpois(1, ticks.observed[1,1])) + 1
data$n.ic <- round(rpois(1, ticks.observed[2,1])) + 1
data$a.ic <- round(rpois(1, ticks.observed[3,1])) + 1

t <- 1
n.adapt <- 2000
n.iter <- 10000
n.chains <- 3
iter2save <- 10000
end <- (length(days)-1)
# end <- 10 # for testing

for (t in 1:end) {
  
  cat("==================================\n")
  cat(round(t/end*100, 2), "% done\n")
  
  ### forecast step ###
  forecast_issue_time <- as.Date(days[t])
  forecast.start.day <- as.character(days[t])  # date forecast issued
  forecast.end.day <- as.character(days[t+1])  # next observation date
  check.day <- as.integer(days[t+1] - days[t]) # day we evaluate forecast and run DA
  
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
  
  data$n.days <- nrow(gdd) # number of days in forecast
  data$gdd <- gdd
  data$met.obs <- obs.temp
  if(any(is.na(obs.temp))){
    data$met.obs.miss <- which(is.na(obs.temp))
    data$met.obs.range <- range(obs.temp, na.rm = TRUE)  
  }
  
  y <- matrix(NA, 3, check.day)
  y[,1] <- ticks.observed[,t]
  data$y <- y
  data$seq.days <- (data$n.days-1):1
  
  # run filter
  enkf <- run_enkf(data, n.adapt = n.adapt, n.chains = 3)
  
  jags.out <- coda.samples(model = enkf$j.model,
                           variable.names = enkf$monitor,
                           n.iter = n.iter)
  
  cat("coda samples done, checking mcmc \n")
  
  ## convergence check
  out <- convergence_check(jags.out, enkf)
  
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
  data <- update_data(parameters, preds) # update with full posterior
  data$l.ic <- round(mean(larva[,ncol(larva)]))
  data$n.ic <- round(mean(nymph[,ncol(nymph)]))
  data$a.ic <- round(mean(adult[,ncol(adult)]))
  
  # thin for saving
  thin <- seq(1, nrow(parameters), length.out = iter2save)
  preds <- preds[thin,]
  parameters <- parameters[thin,]
  
  # save as .RData
  outname <- paste(t, "Tick_hindcast_enkf.RData", sep = "_") # name file
  save(preds, parameters,
       file = file.path(out.dir, outname))
  
}







