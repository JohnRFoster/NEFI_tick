library(mvtnorm)
library(tidyverse)

source("Functions/cary_tick_met_JAGS.R")
source("Functions/site_data_met.R")

# subset data to site.run and met variable of interest
data.all <- cary_ticks_met_JAGS()

all.sites <- c("Green Control", "Henry Control", "Tea Control")
site.dir <- gsub(" Control", "", all.sites)
data <- site_data_met(site = all.sites[1], met.variable = NULL, data.all)

# load in random walk output
load.dir <- paste0("/projectnb/dietzelab/fosterj/FinalOut/A_Correct/RandomWalk/",
                   site.dir[1],
                   "/Combined_thinMat_RandomWalk_",
                   gsub(" ", "", all.sites[1]),
                   ".RData")
load(load.dir)

random_walk <- function(params, N_est, obs, Nmc = 500){
  
  # random samples
  if(Nmc < nrow(params)){
    samps <- sample.int(nrow(params), Nmc)
  } else {
    samps <- 1:nrow(params)
  }
  
  # subset random samples
  params <- params[samps,]
  
  # covariance matrix
  SIGMA <- array(NA, dim = c(3,3,Nmc))
  SIGMA[1,1,] <- params[,"SIGMA[1,1]"]
  SIGMA[1,2,] <- params[,"SIGMA[1,2]"]
  SIGMA[1,3,] <- params[,"SIGMA[1,3]"]
  SIGMA[2,1,] <- params[,"SIGMA[2,1]"]
  SIGMA[2,2,] <- params[,"SIGMA[2,2]"]
  SIGMA[2,3,] <- params[,"SIGMA[2,3]"]
  SIGMA[3,1,] <- params[,"SIGMA[3,1]"]
  SIGMA[3,2,] <- params[,"SIGMA[3,2]"]
  SIGMA[3,3,] <- params[,"SIGMA[3,3]"]
  
  # convert from precision to standard dev
  for(i in 1:Nmc){
    SIGMA[,,i] <- solve(SIGMA[,,i])
  }
  
  pred <- array(NA, dim = c(3, N_est, Nmc))
  for(m in 1:Nmc){
    # pred[,1,m] <- obs[,1]
    if(m %% 500 == 0) cat(m, "mcmc samples complete\n")
    for(t in 1:(N_est-1)){
      # process model
      p <- rmvnorm(1, obs[,t], SIGMA[,,m])
      pred[1, t+1, m] <- max(p[1], 0)
      pred[2, t+1, m] <- max(p[2], 0)
      pred[3, t+1, m] <- max(p[3], 0)
    }  
  }
  return(pred)
}

pred.test <- random_walk(params.mat, data$N_est, data$y)
larva <- t(pred.test[1,,])
nymph <- t(pred.test[2,,])
adult <- t(pred.test[3,,])

quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
larva.q <- apply(larva[,-1], 2, quantile, quants)
nymph.q <- apply(nymph[,-1], 2, quantile, quants)
adult.q <- apply(adult[,-1], 2, quantile, quants)

time <- 1:ncol(larva.q)
plot(time, larva.q[3,], pch = "")
ciEnvelope(time, larva.q[1,], larva.q[5,], col = "lightblue")
ciEnvelope(time, larva.q[2,], larva.q[4,], col = "grey")
lines(time, larva.q[3,])
points(time, data$y[1,-1])


