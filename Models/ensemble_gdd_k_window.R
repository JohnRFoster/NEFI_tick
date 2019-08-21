library(mvtnorm)
library(boot)

source("Functions/Future_Met.R")

ensemble_gdd_k_window <- function(site, params, ic, gdd, N_days, Nmc){
  
  # call data from model run for storage dimensions and indexing
  if(site == "Green Control"){
    s <- 1
  } else if(site == "Henry Control"){
    s <- 2
  } else {
    s <- 3
  }
  
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
  
  SIGMA <- SIGMA.cov <- array(NA, dim = c(3,3,Nmc))
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
    SIGMA.cov[,,i] <- solve(SIGMA[,,i])
  }
  
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
    
    
    
    # zero inflate????
    Nprev <- matrix(NA, 3, 1)
    # Nprev[1,1] <- rpois(1, inflate.larvae)
    # Nprev[2,1] <- rpois(1, inflate.nymph)
    # Nprev[3,1] <- rpois(1, inflate.adult)
    Nprev[1,1] <- ic[1,m]
    Nprev[2,1] <- ic[2,m]
    Nprev[3,1] <- ic[3,m]
    # Nprev[1,1] <- rpois(1, ic[1,m])
    # Nprev[2,1] <- rpois(1, ic[2,m])
    # Nprev[3,1] <- rpois(1, ic[3,m])
    
    ## aggrigate transition matricies
    TRANS <- A[,,1] 
    
    for(t in 1:N_days){
      if(t == 2){
        TRANS <- TRANS %*% A[,,2]
      } else if(t >= 3){
        TRANS <- TRANS %*% A[,,2]
        for(d in 3:t){
          TRANS <- TRANS %*% A[,,d]
        }
      }
      # predict questing ticks
      Ex <- TRANS %*% Nprev
      
      # process error
      est.mvnorm <- rmvnorm(1,Ex,SIGMA.cov[,,m])
      
      # negative truncation and storage
      pred[1,t,m] <- max(est.mvnorm[1], 0)
      pred[2,t,m] <- max(est.mvnorm[2], 0)
      pred[3,t,m] <- max(est.mvnorm[3], 0)
      
      # lar <- max(est.mvnorm[1], 0)
      # nym <- max(est.mvnorm[2], 0)
      # adu <- max(est.mvnorm[3], 0)
      # 
      # pred[1,t,m] <- lar * rbinom(1, 1, params[m,"theta.larvae"])
      # pred[2,t,m] <- nym * rbinom(1, 1, params[m,"theta.nymph"])
      # pred[3,t,m] <- adu * rbinom(1, 1, params[m,"theta.adult"])
      
      Nprev[1,1] <- pred[1,t,m]
      Nprev[2,1] <- pred[2,t,m]
      Nprev[3,1] <- pred[3,t,m]
    }
  }
  return(pred)
} # close function
