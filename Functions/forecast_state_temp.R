forecast_state_temp <- function(type, site, params, ic, N_days, df, Nmc, draw){
  
  # call data from model run for storage dimensions and indexing
  if(site == "Green Control"){
    s <- 1
    t <- 72
  } else if(site == "Henry Control"){
    s <- 2
    t <- 73
  } else {
    s <- 3
    t <- 72
  }

  dt.index <- seq(df, N_days, by = df)
  N_est <- length(dt.index)
  df <- rep(df, length(dt.index))
  
  # get temperature data
  source("Functions/Future_Met.R")
  met.hind <- future_met(site, 10)
  cum.gdd <- met.hind$cum.gdd[1:N_days]
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
  
  alpha.11 <- params[,paste("alpha.11[", s, "]", sep = "")]
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
    alpha.11 <- alpha.11[draw]
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
    alpha.11 <- rep(mean(alpha.11), Nmc)
    alpha.13 <- rep(mean(alpha.13), Nmc)
    alpha.22 <- rep(mean(alpha.22), Nmc)
    alpha.33 <- rep(mean(alpha.33), Nmc)
    alpha.32 <- rep(mean(alpha.32), Nmc)
  }
  
  # select appropraite covariance matrix
  SIGMA <- array(NA, dim = c(3,3,Nmc))
  if("process" %in% type){
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
    SIGMA[1,1,] <- rep(param.mean["SIGMA[1,1]"], Nmc)
    SIGMA[1,2,] <- rep(param.mean["SIGMA[1,2]"], Nmc)
    SIGMA[1,3,] <- rep(param.mean["SIGMA[1,3]"], Nmc)
    SIGMA[2,1,] <- rep(param.mean["SIGMA[2,1]"], Nmc)
    SIGMA[2,2,] <- rep(param.mean["SIGMA[2,2]"], Nmc)
    SIGMA[2,3,] <- rep(param.mean["SIGMA[2,3]"], Nmc)
    SIGMA[3,1,] <- rep(param.mean["SIGMA[3,1]"], Nmc)
    SIGMA[3,2,] <- rep(param.mean["SIGMA[3,2]"], Nmc)
    SIGMA[3,3,] <- rep(param.mean["SIGMA[3,3]"], Nmc)
    for(i in 1:Nmc){
      SIGMA[,,i] <- solve(SIGMA[,,i])
      SIGMA[,,i] <- diag(diag(SIGMA[,,i]))
    }
  }
  
  # random error
  if("random effect" %in% type){
    tau.11 <- 1/sqrt(params[draw, "tau.11"])
    tau.13 <- 1/sqrt(params[draw, "tau.13"])
    tau.22 <- 1/sqrt(params[draw, "tau.22"])
    tau.33 <- 1/sqrt(params[draw, "tau.33"])
    tau.32 <- 1/sqrt(params[draw, "tau.32"])
    alpha.11 <- rnorm(Nmc, 0, tau.11)
    alpha.13 <- rnorm(Nmc, 0, tau.13)
    alpha.22 <- rnorm(Nmc, 0, tau.22)
    alpha.33 <- rnorm(Nmc, 0, tau.33)
    alpha.32 <- rnorm(Nmc, 0, tau.32)
  }
  

  # run mcmc sampling
  for(m in 1:Nmc){
    
    # draw transition matrix
    A[1,1,] <- inv.logit(phi.l.mu[m] + beta.11[m]*met.x + alpha.11[m])
    A[1,3,] <- exp(repro.mu[m] + alpha.13[m])
    A[2,1,] <- inv.logit(grow.ln.mu[m])
    A[2,2,] <- inv.logit(phi.n.mu[m] + beta.22[m]*met.x + alpha.22[m])
    A[3,2,] <- inv.logit(grow.na.mu[m] + alpha.32[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + beta.33[m]*met.x + alpha.33[m])
    
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
      
      ## initial condition
      l <- larva.ic[m]
      n <- nymph.ic[m]
      a <- adult.ic[m]
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
      pred[1,t,m] <- max(est.mvnorm[1], 0)
      pred[2,t,m] <- max(est.mvnorm[2], 0)
      pred[3,t,m] <- max(est.mvnorm[3], 0)
      larva.ic[m] <- pred[1,t,m]
      nymph.ic[m] <- pred[2,t,m]
      adult.ic[m] <- pred[3,t,m]
    }
  }
  return(pred)
} # close function
