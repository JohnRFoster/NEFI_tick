library(mvtnorm)
library(boot)

source("Functions/site_data_met.R")
source("Functions/Future_Met.R")

forecast_state_one_gdd_temp_window <- function(type, site, params, ic, N_est, df, Nmc, draw){
  
  # call data from model run for storage dimensions and indexing
  if(site == "Green Control"){
    s <- 1
    t <- 72
  } else if(site == "Henry Control"){
    s <- 2
    t <- 72
  } else {
    s <- 3
    t <- 72
  }
  
  df <- rep(df, N_est)
  dt.index <- cumsum(df)
  N_days <- dt.index[length(dt.index)]
  
  # get temperature data
  met.hind <- future_met(site, 10)
  gdd <- met.hind$cum.gdd[1:N_days]
  met.x <- met.hind$temp.scale[1:N_days]
  
  # storage
  pred <- array(dim = c(3,N_est,Nmc))
  A <- array(0, dim = c(3,3,N_days))
  
  # mean parameters
  param.mean <- apply(params, 2, mean)
  
  n.prev <- matrix(NA, 3, Nmc)
  
  # select appropriate initial conditions
  if("ic" %in% type){
    larva.ic <- ic[draw,paste("x[1,",t, ",", s,"]",sep="")]
    nymph.ic <- ic[draw,paste("x[2,",t, ",", s,"]",sep="")]
    adult.ic <- ic[draw,paste("x[3,",t, ",", s,"]",sep="")]
  } else {
    IC.mean <- apply(ic, 2, mean)
    larva.ic <- rep(IC.mean[paste("x[1,",t, ",", s,"]",sep="")], Nmc)
    nymph.ic <- rep(IC.mean[paste("x[2,",t, ",", s,"]",sep="")], Nmc)
    adult.ic <- rep(IC.mean[paste("x[3,",t, ",", s,"]",sep="")], Nmc)
  }
  
  IC <- rbind(larva.ic, nymph.ic, adult.ic)
  
  alpha.13 <- params[,paste("alpha.13[", s, "]", sep = "")]
  alpha.22 <- params[,paste("alpha.22[", s, "]", sep = "")]
  alpha.33 <- params[,paste("alpha.33[", s, "]", sep = "")]
  
  # select appropriate parameters
  if("parameter" %in% type){
    phi.l.mu <- params[draw,"phi.l.mu"]
    phi.n.mu <- params[draw,"phi.n.mu"]
    phi.a.mu <- params[draw,"phi.a.mu"]
    grow.ln.mu <- params[draw,"grow.ln.mu"]
    grow.na.mu <- params[draw,"grow.na.mu"]
    repro.mu <- params[draw,"repro.mu"]
    beta.11 <- params[draw, "beta.11"]
    beta.22 <- params[draw, "beta.22"]
    beta.33 <- params[draw, "beta.33"]
    alpha.13 <- alpha.13[draw]
    alpha.22 <- alpha.22[draw]
    alpha.33 <- alpha.33[draw]
  } else {
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
    beta.11 <- rep(param.mean["beta.11"], Nmc)
    beta.22 <- rep(param.mean["beta.22"], Nmc)
    beta.33 <- rep(param.mean["beta.33"], Nmc)
    alpha.13 <- rep(mean(alpha.13), Nmc)
    alpha.22 <- rep(mean(alpha.22), Nmc)
    alpha.33 <- rep(mean(alpha.33), Nmc)
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
  
  # random error
  if("random effect" %in% type){
    tau.13 <- 1/sqrt(params[draw, "tau.13"])
    tau.22 <- 1/sqrt(params[draw, "tau.22"])
    tau.33 <- 1/sqrt(params[draw, "tau.33"])
    alpha.13 <- rnorm(Nmc, 0, tau.13)
    alpha.22 <- rnorm(Nmc, 0, tau.22)
    alpha.33 <- rnorm(Nmc, 0, tau.33)
  }
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    theta.21 <- ifelse((gdd >= 500) & (gdd <= 2000),grow.ln.mu[m],0)
    theta.32 <- ifelse((gdd <= 750) | (gdd >= 2500),grow.na.mu[m],0)
    
    phi.11 <- inv.logit(phi.l.mu[m] + beta.11[m]*met.x)      
    phi.22 <- inv.logit(phi.n.mu[m] + beta.22[m]*met.x + alpha.22[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + beta.33[m]*met.x + alpha.33[m])
    
    A[1,1,] <- phi.11*(1-theta.21) 
    A[2,1,] <- phi.11*theta.21 
    A[2,2,] <- phi.22*(1-theta.32) 
    A[3,2,] <- phi.22*theta.32
    A[1,3,] <- exp(repro.mu[m] + alpha.13[m])
    
    # last observation in time series is first obs for projection
    Nprev <- t(t(IC[,m]))
    
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
} # close function
