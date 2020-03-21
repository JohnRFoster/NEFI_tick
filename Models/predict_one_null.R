#' This is a one step ahead prediction for the GDD_Threshold model
#' 
#' Uncertainty is also partioned 
#' 
#' This model has interacting parameters in the matrix (phi*(1-theta))
#' 
#' Transition is a switch, on when over a threshold or within window, depending
#' on function call
#' 
#' Survival in the questing matrix constant for each life stage      
#'
#' @param type type of uncertainty to simulate. One of "deterministic","ic","parameter","process"
#' @param thresh one of "low" or "window"
#' @param site site, one of "Green Control","Henry Control","Tea Control"
#' @param params matrix of JAGS output of parameters
#' @param ic matrix of JAGS output of states
#' @param data model data fed to JAGS
#' @param Nmc Number of mc samples
#' @param draw Vector of random draws (row numbers) or sampling
#' @export

library(mvtnorm)
library(boot)

source("Functions/site_data_met.R")
source("Functions/ua_parts.R")
source("Functions/obs_prob.R")
source("Functions/build_periodic_matrices_null.R")

predict_one_null <- function(type, site, params, ic, data){
  
  # number of samples
  Nmc <- nrow(params)
  
  data <- site_data_met(site,met.variable=NULL,data)
  
  N_est <- data$N_est
  N_days <- data$N_day
  dt.index <- data$dt.index
  df <- data$df
  gdd <- data$gdd[1:N_days]
  obs.temp <- data$met.obs[dt.index]
  obs.temp[is.na(obs.temp)] <- mean(obs.temp, na.rm = TRUE)
  
  seq.days <- matrix(NA, N_est-1, max(df, na.rm = TRUE))
  for(i in 1:(N_est-1)){
    xx <- (dt.index[i+1]-1):dt.index[i]
    seq.days[i,1:length(xx)] <- xx
  }
  
  # storage
  A <- A.agg <- array(0, dim=c(3,3,N_days))
  TRANS <- array(0, dim=c(3,3,N_est,Nmc))
  pred <- array(dim = c(3,N_est-1,Nmc))
  
  # partition samples
  ua <- ua_parts(params, ic)
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    phi.11 <- inv.logit(ua$phi.l.mu[m])
    phi.22 <- inv.logit(ua$phi.n.mu[m])
    lambda <- ifelse(gdd >= 1500 & gdd <= 2500, ua$repro.mu[m], 0)
    theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, inv.logit(ua$grow.na.mu[m]), 0)
    theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)  
    
    # draw transition matrix
    A[1,1,] <- phi.11*(1-theta.21)
    A[2,1,] <- phi.11*theta.21
    A[2,2,] <- phi.22*(1-theta.32)
    A[3,2,] <- phi.22*theta.32
    A[3,3,] <- inv.logit(ua$phi.a.mu[m])
    A[1,3,] <- lambda
    
    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){                 # loop over the number of sampling days - 1

      for(day in seq.days[t,1]){
        A.agg[1:3,1:3,day] <- A[,,day] %*% A[,,day-1]
      } 
      
      for(day in seq.days[t,2:df[t]]){
        A.agg[1:3,1:3,day] <- A.agg[1:3,1:3,day+1] %*% A[,,day]
      }
      
      ## initial condition
      l <- ua$IC[m, paste("x[1,",t, "]",sep="")]
      n <- ua$IC[m, paste("x[2,",t, "]",sep="")]
      a <- ua$IC[m, paste("x[3,",t, "]",sep="")]
      
      TRANS <- A.agg[,,day]
      
      ## predict new
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      est.mvnorm <- rmvnorm(1,Ex,ua$SIGMA[,,m])
      
      ## Observation model
      
      b.larva <- rbinom(1, 1, ua$theta.larva[m])
      b.nymph <- rbinom(1, 1, ua$theta.nymph[m])
      b.adult <- rbinom(1, 1, ua$theta.adult[m])
      
      pred[1,t,m] <- max(est.mvnorm[1], 0) * b.larva
      pred[2,t,m] <- max(est.mvnorm[2], 0) * b.nymph
      pred[3,t,m] <- max(est.mvnorm[3], 0) * b.adult
      
    }
  }
  return(pred)
}


