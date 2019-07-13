library(mvtnorm)
library(boot)

source("Functions/Future_Met.R")

ensemble_gdd_k_window <- function(site, params, ic, gdd, N_days, Nmc){
  
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
  
  # df <- rep(df, N_est)
  # dt.index <- cumsum(df)
  # N_days <- dt.index[length(dt.index)]
  
  # storage
  pred <- array(dim = c(3,N_days,Nmc))
  A <- array(0, dim = c(3,3,N_days))

  alpha.13 <- params[, paste("alpha.13[", s, "]", sep = "")]
  alpha.11 <- params[, paste("alpha.11[", s, "]", sep = "")]
  alpha.21 <- params[, paste("alpha.21[", s, "]", sep = "")]
  alpha.33 <- params[, paste("alpha.33[", s, "]", sep = "")]
  phi.l.mu <- params[,"phi.l.mu"]
  phi.n.mu <- params[,"phi.n.mu"]
  phi.a.mu <- params[,"phi.a.mu"]
  grow.ln.mu <- params[,"grow.ln.mu"]
  grow.na.mu <- params[,"grow.na.mu"]
  repro.mu <- params[,"repro.mu"]
  k.l2n.low <- params[, "k.l2n.low"]
  k.l2n.high <- params[, "k.l2n.high"]
  k.n2a.low <- params[, "k.n2a.low"]
  k.n2a.high <- params[, "k.n2a.high"]
  theta.larvae <- params[, "theta.larvae"]
  theta.nymph <- params[, "theta.nymph"]
  theta.adult <- params[, "theta.adult"]

  SIGMA <- array(NA, dim = c(3,3,Nmc))
  SIGMA[1,1,] <- params[, "SIGMA[1,1]"]
  SIGMA[1,2,] <- params[, "SIGMA[1,2]"]
  SIGMA[1,3,] <- params[, "SIGMA[1,3]"]
  SIGMA[2,1,] <- params[, "SIGMA[2,1]"]
  SIGMA[2,2,] <- params[, "SIGMA[2,2]"]
  SIGMA[2,3,] <- params[, "SIGMA[2,3]"]
  SIGMA[3,1,] <- params[, "SIGMA[3,1]"]
  SIGMA[3,2,] <- params[, "SIGMA[3,2]"]
  SIGMA[3,3,] <- params[, "SIGMA[3,3]"]
  
  # convert from precision to standard dev
  for(i in 1:Nmc){
    SIGMA[,,i] <- solve(SIGMA[,,i])
  }
  
  # random error
  # tau.13 <- 1/sqrt(params[draw, "tau.13"])
  # tau.21 <- 1/sqrt(params[draw, "tau.21"])
  # tau.33 <- 1/sqrt(params[draw, "tau.33"])
  # alpha.13 <- rnorm(Nmc, 0, tau.13)
  # alpha.21 <- rnorm(Nmc, 0, tau.21)
  # alpha.33 <- rnorm(Nmc, 0, tau.33)
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    t21 <- inv.logit(grow.ln.mu[m] + alpha.21[m])
    
    theta.21 <- ifelse((gdd >= k.l2n.low[m]) & (gdd <= k.l2n.high[m]),t21,0)
    theta.32 <- ifelse((gdd <= k.n2a.low[m]) | (gdd >= k.n2a.high[m]),grow.na.mu[m],0)
    
    phi.11 <- inv.logit(phi.l.mu[m])
    phi.22 <- inv.logit(phi.n.mu[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + alpha.33[m])
    
    A[1,1,] <- phi.11*(1-theta.21) 
    A[2,1,] <- phi.11*theta.21 
    A[2,2,] <- phi.22*(1-theta.32) 
    A[3,2,] <- phi.22*theta.32
    A[1,3,] <- exp(repro.mu[m] + alpha.13[m])
    
    ## binary outcome of observation by life stage and
    ## blend the poisson and zero inflation models
    inflate.larvae <- ic[1,1]*rbinom(1, 1, theta.larvae[m]) + 0.1
    inflate.nymph <- ic[2,1]*rbinom(1, 1, theta.nymph[m]) + 0.1
    inflate.adult <- ic[3,1]*rbinom(1, 1, theta.adult[m]) + 0.1
    
    # zero inflate????
    Nprev <- matrix(NA, 3, 1)
    Nprev[1,1] <- rpois(1, inflate.larvae)
    Nprev[2,1] <- rpois(1, inflate.nymph)
    Nprev[3,1] <- rpois(1, inflate.adult)
    # Nprev[1,1] <- rpois(1, ic[1,1])
    # Nprev[2,1] <- rpois(1, ic[2,1])
    # Nprev[3,1] <- rpois(1, ic[3,1])
    
    ## aggrigate transition matricies
    TRANS <- A[,,1] %*% A[,,2]
    
    for(t in 3:N_days){
      TRANS <- TRANS %*% A[,,d]
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
  return(pred)
} # close function
