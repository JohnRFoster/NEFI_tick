library(mvtnorm)
library(tidyverse)
library(ecoforecastR)
library(lubridate)
library(gridExtra)
library(scoringRules)

source("Functions/cary_tick_met_JAGS.R")
source("Functions/site_data_met.R")
source("Functions/model_performance_metrics.R")


# subset data to site.run and met variable of interest
data.all <- cary_ticks_met_JAGS()

all.sites <- c("Green Control", "Henry Control", "Tea Control")
site.dir <- gsub(" Control", "", all.sites)
data <- site_data_met(site = all.sites[1], met.variable = NULL, data.all)
dates <- read.csv("/projectnb/dietzelab/fosterj/Data/tick_cleaned") %>% 
  filter(Grid == all.sites[1]) %>%
  pull(DATE)
  
# load in random walk output
load.dir <- paste0("/projectnb/dietzelab/fosterj/FinalOut/A_Correct/RandomWalk/",
                   site.dir[1],
                   "/Combined_thinMat_RandomWalk_",
                   gsub(" ", "", all.sites[1]),
                   ".RData")
load(load.dir)

random_walk <- function(params, N_est, start, obs, Nmc = 500){
  
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
  
  pred.seq <- start:(N_est-1)
  pred <- array(NA, dim = c(3, N_est, Nmc))
  for(m in 1:Nmc){
    pred[,start,m] <- obs[,start]
    if(m %% 500 == 0) cat(m, "mcmc samples complete\n")
    for(t in pred.seq){
      # process model
      p <- rmvnorm(1, pred[,t,m], SIGMA[,,m])
      pred[1, t+1, m] <- max(p[1], 0)
      pred[2, t+1, m] <- max(p[2], 0)
      pred[3, t+1, m] <- max(p[3], 0)
    }  
  }
  return(pred)
}

k.end <- 50
pred.k <- larva.q <- nymph.q <- adult.q <- list()
crps.l <- crps.n <- crps.a <- matrix(NA, k.end, data$N_est)
quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
for(k in 1:k.end){
  print(k)
  pred.k[[k]] <- random_walk(params.mat, data$N_est, k, data$y, Nmc = 100)
  
  remove.seq <- 1:k
  
  # continuous ranked probability score
  crps.l[k, -remove.seq] <- crps_sample(data$y[1,-remove.seq], pred.k[[k]][1,-remove.seq,])
  crps.n[k, -remove.seq] <- crps_sample(data$y[2,-remove.seq], pred.k[[k]][2,-remove.seq,])
  crps.a[k, -remove.seq] <- crps_sample(data$y[3,-remove.seq], pred.k[[k]][3,-remove.seq,])

  larva <- t(pred.k[[k]][1,,])
  nymph <- t(pred.k[[k]][2,,])
  adult <- t(pred.k[[k]][3,,])

  larva.q[[k]] <- apply(larva[,-remove.seq], 2, quantile, quants)
  nymph.q[[k]] <- apply(nymph[,-remove.seq], 2, quantile, quants)
  adult.q[[k]] <- apply(adult[,-remove.seq], 2, quantile, quants)
}


crps.l.mu <- apply(crps.l, 2, mean, na.rm = TRUE)
crps.n.mu <- apply(crps.n, 2, mean, na.rm = TRUE)
crps.a.mu <- apply(crps.a, 2, mean, na.rm = TRUE)

plot(crps.l.mu)
plot(crps.n.mu)
plot(crps.a.mu)

plot_time_series <- function(quants, obs, k){
  col <- c('#a6bddb', '#2b8cbe')
  alpha <- 1
  axis.size <- 18
  axis.title <- 20
  
  quants %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(date = ymd(dates[(k+1):length(dates)]),
           obs = obs[1,(k+1):length(dates)]) %>% 
    ggplot(aes(x = date, group = 1)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%"), alpha = alpha) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = "25% - 75%"), alpha = alpha) +
    geom_line(aes(y = `50%`, fill = "Median"), size = 1) +
    geom_point(aes(y = obs), size = 2) +
    scale_fill_manual(name = "Prediction Interval",
                      values = c("2.5% - 97.5%" = col[1],
                                 "25% - 75%" = col[2],
                                 "Median" = "black"),
                      breaks = c("2.5% - 97.5%",
                                 "25% - 75%",
                                 "Median")) +
    theme_classic() +
    theme(axis.text = element_text(size = axis.size),
          axis.title = element_text(size = axis.title),
          plot.title = element_text(size = axis.title),
          legend.position = "none")
}

  
  


   

par(mfrow = c(1,3))
for(k in 1:5){
  l <- plot_time_series(larva.q[[k]], data$y, k) +
    labs(title = "Larvae",
         x = "",
         y = "Individuals") +
    coord_cartesian(ylim = c(0, 6000))
  
  n <- plot_time_series(nymph.q[[k]], data$y, k) +
    labs(title = "Nymph",
         x = "",
         y = "Individuals") +
    coord_cartesian(ylim = c(0, 200))
  
  a <- plot_time_series(adult.q[[k]], data$y, k) +
    labs(title = "Adult",
         x = "",
         y = "Individuals") +
    coord_cartesian(ylim = c(0, 100)) 
  
  grid.arrange(l, n, a, nrow = 1)
}





str(larva.q[[5]])






