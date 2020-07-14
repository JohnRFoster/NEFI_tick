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



predict_one <- function(type, site, met.variable, params, ic, data, process.type,
                        A_func = NULL, A = NULL, dir = "", HB = FALSE){
  
  source(paste0(dir, "Functions/site_data_met.R"))
  # source(paste0(dir, "Functions/cary_tick_met_JAGS.R"))
  source(paste0(dir, "Functions/ua_parts.R"))
  source(paste0(dir, "Functions/obs_prob.R"))
  source(paste0(dir, "Functions/get_month_effect.R"))
  
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
  ua <- ua_parts(params, ic, type, process.type)
  
  # get A
  if(is.null(A)){
    A <- A_func(ua, data)  
  }
  
  obs.prob <- obs_prob(ua, N_est, obs.temp)
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){    # loop over the number of sampling days - 1
      
      for(day in seq.days[t,1]){
        A.agg[1:3,1:3,day] <- A[,,day,m] %*% A[,,day-1,m]
      } 
      
      for(day in seq.days[t,2:df[t]]){
        A.agg[1:3,1:3,day] <- A.agg[1:3,1:3,day+1] %*% A[,,day,m]
      }
 
      ## initial condition
      if(is.numeric(HB)){
        end <- paste0(",", HB, "]")
        l <- ua$IC[m, paste0("x[1,",t,end)]
        n <- ua$IC[m, paste0("x[2,",t,end)]
        a <- ua$IC[m, paste0("x[3,",t,end)]  
      } else {
        l <- ua$IC[m, paste0("x[1,",t,"]")]
        n <- ua$IC[m, paste0("x[2,",t,"]")]
        a <- ua$IC[m, paste0("x[3,",t,"]")]  
      }
      
      TRANS <- A.agg[,,day]
      
      ## predict new
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      
      # check for month effect in jags output and add to Ex
      if(any(grepl("alpha.month", names(ua)))){
        month.effect <- get_month_effect(t, m, ua, data$month.index)
        Ex <- Ex + month.effect
      }
      
      est.mvnorm <- rep(NA, 3)
      if(process.type == "multi"){
        est.mvnorm <- rmvnorm(1,Ex,ua$SIGMA[,,m])  
      } else if(process.type == "independent"){
        est.mvnorm[1] <- rnorm(1,Ex[1],ua$tau.l[m])  
        est.mvnorm[2] <- rnorm(1,Ex[2],ua$tau.n[m])  
        est.mvnorm[3] <- rnorm(1,Ex[3],ua$tau.a[m])  
      }
      
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

# load("../FinalOut/A_Correct/NULL/Green/Combined_thinMat_NULL_GreenControl.RData")
# source("Functions/build_periodic_matrices_null.R")
# data <- cary_ticks_met_JAGS()
# site <- "Green Control"
# type <- c("ic", "process", "parameter")
# Nmc <- 300
# draw <- sample.int(nrow(params.mat), Nmc, replace = TRUE)
# params <- params.mat[draw,]
# ic <- predict.mat[draw,]
# 
# test <- predict_one("parameter",
#                     site,
#                     NULL,
#                     params,
#                     ic,
#                     data,
#                     build_periodic_matrices_null)
# 
# 
# load("../FinalOut/Independent_Fits/GDDThreshold/RH_ObsProc/beta_111/Green/Combined_thinMat_MaxRH_ObsProc_beta_111_K_set_GreenControl.RData")
# source("Functions/build_periodic_matrices_1proc.R")
# data <- cary_ticks_met_JAGS()
# site <- "Green Control"
# type <- c("ic", "process", "parameter")
# Nmc <- 300
# draw <- sample.int(nrow(params.mat), Nmc, replace = TRUE)
# params <- params.mat[draw,]
# ic <- predict.mat[draw,]
# test <- predict_one(type,
#                     site,
#                     "max rh",
#                     params,
#                     ic,
#                     data,
#                     build_periodic_matrices_1proc)
