#' This is the forecast function for the GDD_Threshold (null) model
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
#' @param N_est number of future (sampling) occasians to forecast
#' @param df the number of days between future (sampling) occasions
#' @param Nmc Number of mc samples
#' @param draw Vector of random draws (row numbers) or sampling
#' @export

library(mvtnorm)

# source("Functions/site_data_met.R")
# source("Functions/Future_Met.R")
# source("Functions/get_last_day.R")

forecast_state_gdd_null <- function(type, thresh, site, params, ic, N_est, df, Nmc, draw){
  
  # site index
  if(site == "Green Control"){
    s <- 1
  } else if(site == "Henry Control"){
    s <- 2
  } else {
    s <- 3
  }
  
  df <- rep(df, N_est)
  dt.index <- cumsum(df)
  N_days <- dt.index[length(dt.index)]
  
  # future met - technically hindcasting
  met.hind <- future_met(site, 10)
  gdd <- met.hind$cum.gdd[1:N_days]
  
  # select appropriate initial conditions
  last <- get_last_day()
  n.obs <- last[s,"n.obs"]
  ic <- ic[,grep(paste(n.obs, "]", sep = ""), colnames(ic))]

  # storage
  pred <- array(dim = c(3,N_est,Nmc))
  n.prev <- matrix(NA, 3, Nmc)
  A <- array(0, dim=c(3,3,N_days))
  
  # mean parameters
  param.mean <- apply(params, 2, mean)
  
  # select appropriate initial conditions
  if("ic" %in% type){
    IC <- ic[draw,]
  } else {
    IC.mean <- apply(ic, 2, mean)
    IC.names <- names(IC.mean)
    IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
    colnames(IC) <- IC.names
  }
  
  # select appropriate parameters
  if("parameter" %in% type){
    phi.l.mu <- params[draw,"phi.l.mu"]
    phi.n.mu <- params[draw,"phi.n.mu"]
    phi.a.mu <- params[draw,"phi.a.mu"]
    grow.ln.mu <- params[draw,"grow.ln.mu"]
    grow.na.mu <- params[draw,"grow.na.mu"]
    repro.mu <- params[draw,"repro.mu"]
  } else {
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
  }
  
  # select appropraite covariance matrix
  if("process" %in% type){
    SIGMA <- array(NA, dim = c(3,3,Nmc))
    SIGMA[1,1,] <- params[draw,"SIGMA[1,1]"]
    SIGMA[1,2,] <- params[draw,"SIGMA[1,2]"]
    SIGMA[1,3,] <- params[draw,"SIGMA[1,3]"]
    SIGMA[2,1,] <- params[draw,"SIGMA[2,1]"]
    SIGMA[2,2,] <- params[draw,"SIGMA[2,2]"]
    SIGMA[2,3,] <- params[draw,"SIGMA[2,3]"]
    SIGMA[3,1,] <- params[draw,"SIGMA[3,1]"]
    SIGMA[3,2,] <- params[draw,"SIGMA[3,2]"]
    SIGMA[3,3,] <- params[draw,"SIGMA[3,3]"]
    
    # convert from precision to standard dev
    for(i in 1:Nmc){
      SIGMA[,,i] <- solve(SIGMA[,,i])
    }
  } else {
    SIGMA <- array(0, dim = c(3,3,Nmc))
  }
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    phi.11 <- inv.logit(phi.l.mu[m])
    phi.22 <- inv.logit(phi.n.mu[m])
    theta.32 <- ifelse((gdd <= 750) | (gdd >= 2500),grow.na.mu[m],0)
    
    if(thresh == "low"){
      theta.21 <- ifelse((gdd >= 500),grow.ln.mu[m],0)   
    } else {
      theta.21 <- ifelse((gdd >= 500) & (gdd <= 2000),grow.ln.mu[m],0)  
    }
    
    # draw transition matrix
    A[1,1,] <- phi.11*(1-theta.21)
    A[2,1,] <- phi.11*theta.21
    A[2,2,] <- phi.22*(1-theta.32)
    A[3,2,] <- phi.22*theta.32
    A[3,3,] <- inv.logit(phi.a.mu[m])
    A[1,3,] <- exp(repro.mu[m])
    
    # last observation in time series is first obs for projection
    Nprev <- t(t(IC[m,]))         
    
    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){                 # loop over the number of sampling days - 1
      if(t == 1){
        TRANS <- A[,,1] %*% A[,,2]
        for(d in 3:df[1]){
          TRANS <- TRANS %*% A[,,d]
        }
      } else {
        for(d in 1:df[t]){
          if(d == 1){
            TRANS <- A[,,dt.index[t-1]] %*% A[,,dt.index[t-1]+1]
          } else {
            TRANS <- TRANS %*% A[,,dt.index[t-1]+d]
          }
        }
      }
      
      # predict questing ticks
      Ex <- TRANS %*% Nprev
      
      # process error
      est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
      
      # negative truncation and storage
      pred[1,t,m] <- max(est.mvnorm[1], 0)
      pred[2,t,m] <- max(est.mvnorm[2], 0)
      pred[3,t,m] <- max(est.mvnorm[3], 0)
      
      Nprev[1,1] <- pred[1,t,m]
      Nprev[2,1] <- pred[2,t,m]
      Nprev[3,1] <- pred[3,t,m]
    }
  }
  return(pred)
}