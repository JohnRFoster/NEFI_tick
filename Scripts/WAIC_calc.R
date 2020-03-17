library(ecoforecastR)

source("Functions/site_data_met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/ua_parts.R")
source("Functions/obs_prob.R")

data.hb <- cary_ticks_met_JAGS()

site.folders <- c("Green", "Henry", "Tea")
sites <- c("Green Control", "Henry Control", "Tea Control")
top.dir <- "../FinalOut/A_Correct/ObsModel/L1.N1.A1/"
met.variable <- NULL

for(i in 1:3){
  dir <- paste0(top.dir, site.folders[i])
  model <- paste0("Combined_thinMat_Obs_111_Beta_", site.folders[i], "Control.RData")
  load(file.path(dir, model))
  
  cat("Calculating WAIC for", model, "\n")

  data <- site_data_met(site = sites[i], 
                        met.variable = met.variable,
                        data = data.hb)
  obs.temp <- data$met.obs
  obs.temp[data$met.obs.miss] <- mean(obs.temp, na.rm = TRUE)
  N_est <- data$N_est
  y <- data$y
  ua <- ua_parts(params.mat, predict.mat, c("parameter", "ic"))
  
  # storage
  like.l <- like.n <- like.a <- matrix(NA, nrow(predict.mat), N_est)
  
  for(t in 1:N_est){ # loop over observations
    
    # observation probability
    obs.prob <- obs_prob(ua, obs.temp[t])
    
    # zero inflate latent state
    num <- length(obs.prob$theta.larva)
    l <- ua$IC[, paste("x[1,",t, "]",sep="")] * rbinom(num, 1, obs.prob$theta.larva)
    n <- ua$IC[, paste("x[2,",t, "]",sep="")] * rbinom(num, 1, obs.prob$theta.nymph)
    a <- ua$IC[, paste("x[3,",t, "]",sep="")] * rbinom(num, 1, obs.prob$theta.adult)
    
    # calculate likelihood
    like.l[,t] <- dpois(y[1,t], l)
    like.n[,t] <- dpois(y[2,t], n)
    like.a[,t] <- dpois(y[3,t], a)
  }
  
  like <- cbind(like.l, like.n, like.n)
  fbar <- colMeans(like)
  like[like==0] <- 1E-10
  Pw <- sum(apply(log(like),2,var))
  WAIC <- -2*sum(log(fbar))+2*Pw
  
  cat("Pw:", Pw, "\n")
  cat("WAIC:", WAIC, "\n\n\n")
}








