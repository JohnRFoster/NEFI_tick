source("Functions/Future_Met.R")
source("Functions/get_last_day.R")

forecast_with_molting_precip <- function(type, site, params, ic, N_est, df, threshold, Nmc, draw){
  
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
  
  last <- get_last_day()
  
  n.obs <- last[s,"n.obs"]
  
  met.hind <- future_met(site, 10)
  cum.gdd <- met.hind$cum.gdd[1:N_days]
  precip <- met.hind$precip[1:N_days]

  # storage
  pred <- array(dim = c(3,N_est,Nmc))
  A <- array(0, dim = c(3,3,N_days))
  
  # mean parameters
  param.mean <- apply(params, 2, mean)
  
  # select appropriate initial conditions
  ic <- ic[,grep(paste(n.obs, ",", s, "]", sep = ""), colnames(ic))]

  if("ic" %in% type){
    IC <- ic[draw,]
  } else {
    IC.mean <- apply(ic, 2, mean)
    IC.names <- names(IC.mean)
    IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
    colnames(IC) <- IC.names
  }
  
  alpha.13 <- params[,paste("alpha.13[", s, "]", sep = "")]
  alpha.22 <- params[,paste("alpha.22[", s, "]", sep = "")]
  alpha.33 <- params[,paste("alpha.33[", s, "]", sep = "")]
  alpha.32 <- params[,paste("alpha.32[", s, "]", sep = "")]
  
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
    alpha.32 <- alpha.32[draw]
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
    alpha.32 <- rep(mean(alpha.32), Nmc)
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
  } else if("process mean" %in% type){
    SIGMA <- array(NA, dim = c(3,3,Nmc))
    SIGMA[1,1,] <- rep(param.mean["SIGMA[1,1]"], Nmc)
    SIGMA[1,2,] <- rep(param.mean["SIGMA[1,2]"], Nmc)
    SIGMA[1,3,] <- rep(param.mean["SIGMA[1,3]"], Nmc)
    SIGMA[2,1,] <- rep(param.mean["SIGMA[2,1]"], Nmc)
    SIGMA[2,2,] <- rep(param.mean["SIGMA[2,2]"], Nmc)
    SIGMA[2,3,] <- rep(param.mean["SIGMA[2,3]"], Nmc)
    SIGMA[3,1,] <- rep(param.mean["SIGMA[3,1]"], Nmc)
    SIGMA[3,2,] <- rep(param.mean["SIGMA[3,2]"], Nmc)
    SIGMA[3,3,] <- rep(param.mean["SIGMA[3,3]"], Nmc)
    
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
    tau.32 <- 1/sqrt(params[draw, "tau.32"])
    alpha.13 <- rnorm(Nmc, 0, tau.13)
    alpha.22 <- rnorm(Nmc, 0, tau.22)
    alpha.33 <- rnorm(Nmc, 0, tau.33)
    alpha.32 <- rnorm(Nmc, 0, tau.32)
  }
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    # draw transition matrix
    A[1,1,] <- inv.logit(phi.l.mu[m] + beta.11[m]*precip)
    A[1,3,] <- exp(repro.mu[m] + alpha.13[m])
    A[2,1,] <- inv.logit(grow.ln.mu[m])
    A[2,2,] <- inv.logit(phi.n.mu[m] + beta.22[m]*precip + alpha.22[m])
    A[3,2,] <- inv.logit(grow.na.mu[m] + alpha.32[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + beta.33[m]*precip + alpha.33[m])
    
    ## aggrigate transition matricies
    larva.molt <- nymph.molt <- adult.molt <- rep(0, length = 100) # initialize
    Nprev <- t(t(IC[m,]))         # last observation in time series is first obs for projection
    
    for(t in 1:(N_est)){   # loop over the number of sampling days
      if(t == 1){
        TRANS <- A[,,1] %*% A[,,2]
        for(d in 3:(df[1]-2)){
          TRANS <- TRANS %*% A[,,d]
        }
      } else {
        for(d in 1:df[t]){
          if(d == 1){
            TRANS <- A[,,dt.index[t-1]] %*% A[,,dt.index[t-1]+1]
          } else {
            TRANS <- TRANS %*% A[,,dt.index[t-1]+d-1]
          }
        }
      }
      
      # predict questing ticks
      Ex <- TRANS %*% Nprev
      # est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
      # pred[1,t,m] <- max(est.mvnorm[1], 0)
      # pred[2,t,m] <- max(est.mvnorm[2], 0)
      # pred[3,t,m] <- max(est.mvnorm[3], 0)
      
      # number transitioning
      nymph.molt[t] <- Nprev[1,1] * TRANS[2,1]
      adult.molt[t] <- Nprev[2,1] * TRANS[3,2]
      larva.molt[t] <- Nprev[3,1] * TRANS[1,3]

      # add transition numbers to questing numbers if GDD is reached
      # new larvae
      if(cum.gdd[dt.index[t]] >= threshold[1] & cum.gdd[dt.index[t]] <= 2350){
        # pred[1,t,m] <- sum(larva.molt, pred[1,t,m])
        Ex[1,1] <- sum(larva.molt, Ex[1,1])
        larva.molt <- rep(0, length = 100) # reset
      }

      # new nymphs
      if(cum.gdd[dt.index[t]] >= threshold[2] & cum.gdd[dt.index[t]] <= 2000){
        # pred[2,t,m] <- sum(nymph.molt, pred[2,t,m])
        Ex[2,1] <- sum(nymph.molt, Ex[2,1])
        nymph.molt <- rep(0, length = 100) # reset
      }

      # new adults
      if(!(cum.gdd[dt.index[t]] >= 750 & cum.gdd[dt.index[t]] <= 2500)){
        # pred[3,t,m] <- sum(adult.molt, pred[3,t,m])
        Ex[3,1] <- sum(adult.molt, Ex[3,1])
        adult.molt <- rep(0, length = 100) # reset
      }
      
      est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
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
