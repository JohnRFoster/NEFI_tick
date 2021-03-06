
tick_forecast_parallel <- function(Nmc.seq){
  
  library(boot)
  library(mvtnorm)
  
  # number of samples
  Nmc <- length(Nmc.seq)
  
  params <- params[Nmc.seq,]
  
  # storage
  pred <- array(dim = c(3,n.days+1,Nmc))
  pred[,1,] <- t(ic[Nmc.seq,])
  
  # partition samples
  ua <- ua_parts(params, ic, c("process", "parameter"))
  obs.prob <- obs_prob(ua, n.days, obs.temp[,1])
  
  # build A  matrices
  A <- array(0, dim=c(3,3,n.days)) # storage
  for(m in seq_along(ua$phi.a.mu)){
    
    # daily survival for each life stage, given met driver or not
    if("beta.l" %in% names(ua)){
      phi.11 <- boot::inv.logit(ua$phi.l.mu[m] + ua$beta.l[m]*met.data[,1])
    } else {
      phi.11 <- boot::inv.logit(ua$phi.l.mu[m])
    }
    
    if("beta.n" %in% names(ua)){
      phi.22 <- boot::inv.logit(ua$phi.n.mu[m] + ua$beta.n[m]*met.data[,1])
    } else {
      phi.22 <- boot::inv.logit(ua$phi.n.mu[m])
    }
    
    if("beta.a" %in% names(ua)){
      phi.33 <- boot::inv.logit(ua$phi.a.mu[m] + ua$beta.a[m]*met.data[,1])
    } else {
      phi.33 <- boot::inv.logit(ua$phi.a.mu[m])
    }
    
    if("rho.l" %in% names(ua)){
      larva.start <- ua$rho.l
    } else {
      larva.start <- rep(1500, Nmc)
    }
    
    if("rho.n" %in% names(ua)){
      nymph.start <- ua$rho.n
    } else {
      nymph.start <- rep(500, Nmc)
    }
    
    if("rho.a" %in% names(ua)){
      adult.start <- ua$rho.n
    } else {
      adult.start <- rep(2500, Nmc)
    }
    
    # reproduction and transition
    lambda <- ifelse(gdd[,1] >= larva.start[m] & gdd[,1] <= 2500, ua$repro.mu[m], 0) 
    theta.32 <- ifelse(gdd[,1] <= 1000 | gdd[,1] >= adult.start[m], inv.logit(ua$grow.na.mu[m]), 0) 
    theta.21 <- ifelse(gdd[,1] >= nymph.start[m] & gdd[,1] <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)  
    
    # build transition matrix
    A[1,1,] <- phi.11*(1-theta.21)
    A[2,1,] <- phi.11*theta.21
    A[2,2,] <- phi.22*(1-theta.32)
    A[3,2,] <- phi.22*theta.32
    A[3,3,] <- phi.33
    A[1,3,] <- lambda
    
    
    for(t in 1:n.days){    
      
      # aggregate transition matrices
      if (t == 1){
        TRANS <- A[,,1]
      } else if (t == 2){
        TRANS <- A[,,2] %*% A[,,1]
      } else {
        TRANS <- A[,,t]
        for(day in t:2){
          TRANS <- TRANS %*% A[,,day-1]
        }
      }
      
      # initial condition
      l <- ic[m, 1]
      n <- ic[m, 2]
      a <- ic[m, 3]
      
      # predict new
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      est.mvnorm <- rmvnorm(1,Ex,ua$SIGMA[,,m])
      
      # Observation model
      pred[1,t+1,m] <- max(est.mvnorm[1], 0)
      pred[2,t+1,m] <- max(est.mvnorm[2], 0)
      pred[3,t+1,m] <- max(est.mvnorm[3], 0)
    }
  }
  pred <- round(pred) # ticks are discrete
  return(pred)
}

