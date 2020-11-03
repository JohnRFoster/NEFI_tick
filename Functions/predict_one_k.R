#' This is a one step ahead prediction for tick models
#' 
#' Uncertainty is also partioned 
#' 
#' This model has interacting parameters in the matrix (phi*(1-theta))
#' 
#' Transition is a switch, on when over a threshold or within window, depending
#' on function call
#' 
#'
#' @param type type of uncertainty to simulate. One of "deterministic","ic","parameter","process"
#' @param site site, one of "Green Control","Henry Control","Tea Control"
#' @param met.variable met driver in process model
#' @param params matrix of JAGS output of parameters
#' @param ic matrix of JAGS output of states
#' @param data model data fed to JAGS
#' @param A_func function that specifies how to build A matrices, default = NULL
#' @param A A matrix, default = NULL. Must specify either A_func OR A
#' @export

library(mvtnorm)
library(boot)



predict_one_k <- function(start, type, site, met.variable, params, ic, data, 
                          use.mice = FALSE, dir = "", HB = FALSE){
  
  source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/site_data_met.R")
  source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/ua_parts.R")
  source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/obs_prob.R")
  source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/get_month_effect.R")
  
  # number of samples
  Nmc <- nrow(params)
  
  data <- site_data_met(site, met.variable, data, dir = dir, time.effect = "month")
  
  N_est <- data$N_est
  N_days <- data$N_day
  dt.index <- data$dt.index
  df <- data$df
  gdd <- data$gdd[1:N_days]
  obs.temp <- data$met.obs
  obs.temp[is.na(obs.temp)] <- mean(obs.temp, na.rm = TRUE)
  obs.temp <- obs.temp[1:N_est]
  
  seq.days <- matrix(NA, N_est-1, max(df, na.rm = TRUE))
  for(i in 1:(N_est-1)){
    xx <- (dt.index[i+1]-1):dt.index[i]
    seq.days[i,1:length(xx)] <- xx
  }
  
  # storage
  A.agg <- array(0, dim=c(3,3,N_days))
  pred <- array(dim = c(3,N_est-1,Nmc))
  
  # partition samples
  ua <- ua_parts(params, ic, type)
  
  if(use.mice | !is.null(met.variable)){
    cat("build_periodic_matrices_1proc\n")
    source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/build_periodic_matrices_1proc.R")
    mice.site <- ifelse(use.mice, site, NA)
    A <- build_periodic_matrices_1proc(ua, data, mice.site = mice.site)
  } else {
    cat("build_periodic_matrices_null\n")
    source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/build_periodic_matrices_null.R")
    A <- build_periodic_matrices_null(ua, data, HB)
  }
  
  obs.prob <- obs_prob(ua, N_est, obs.temp)
  
  # run mcmc sampling
  cat("Starting mcmc sampling\n")
  for(m in 1:Nmc){
    
    ## aggrigate transition matricies
    for(t in start:(N_est-1)){    # loop over the number of sampling days - 1
      
      for(day in seq.days[t,1]){
        A.agg[1:3,1:3,day] <- A[,,day,m] %*% A[,,day-1,m]
      } 
      
      for(day in seq.days[t,2:df[t]]){
        A.agg[1:3,1:3,day] <- A.agg[1:3,1:3,day+1] %*% A[,,day,m]
      }
      
      ## initial condition
      if(t == start){
        if(is.numeric(HB)){
          end <- paste0(",", HB, "]")
          l <- ua$IC[m, paste0("x[1,",start,end)]
          n <- ua$IC[m, paste0("x[2,",start,end)]
          a <- ua$IC[m, paste0("x[3,",start,end)]  
        } else {
          l <- ua$IC[m, paste0("x[1,",start,"]")]
          n <- ua$IC[m, paste0("x[2,",start,"]")]
          a <- ua$IC[m, paste0("x[3,",start,"]")]  
        }  
        obs <- as.matrix(c(l,n,a),3,1)
      } else {
        obs <- pred[,t-1,m]
      }
      
      
      TRANS <- A.agg[,,day]
      
      ## predict new
      Ex <- TRANS %*% obs
      
      # check for month effect in jags output and add to Ex
      if(any(grepl("alpha.month", names(ua)))){
        month.effect <- get_month_effect(t, m, ua, data$month.index)
        Ex <- Ex + month.effect
      }
      
      est.mvnorm <- rep(NA, 3)
      est.mvnorm <- rmvnp(1,as.vector(Ex),ua$SIGMA[,,m])  
      
      b.larva <- rbinom(1, 1, obs.prob$theta.larva[m,t+1])
      b.nymph <- rbinom(1, 1, obs.prob$theta.nymph[m,t+1])
      b.adult <- rbinom(1, 1, obs.prob$theta.adult[m,t+1])
      
      ## Observation model
      pred[1,t,m] <- max(est.mvnorm[1], 0)*b.larva #+ 1E-10 
      pred[2,t,m] <- max(est.mvnorm[2], 0)*b.nymph #+ 1E-10 
      pred[3,t,m] <- max(est.mvnorm[3], 0)*b.adult #+ 1E-10 
      
    }
  }
  return(pred)
}

