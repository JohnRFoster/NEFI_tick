##   This is the HB model for Cary Ticks on a Daily scale.  
##
##   Random site effects on adult and larvae survival, larvae-to-nymph transition,
##   reproduction, and the upper and lower larvae-to-nymph transtion threshold
##
##   This model has interacting parameters in the matrix (phi*(1-psi))
##   
##   Transition is a threshold window
##
##   Survival in the questing matrix constant for each life stage        
##                                                                 
##   The state is estimated for every sampling day only, 
##   but demographic params are estimated daily.                       
##
##   Zero Inflated Poisson Data model by life stage. 

#'@param n.adapt number of adaption iterations
#'@param n.chains number of chains, 1 for parallel batch jobs
#'@param burnin number of burnin iterations
#'@param thin thinning interval
#'@param n.iter number of iterations post burnin.     

library(ecoforecastR)

source("Functions/cary_tick_met_JAGS.R")
run_model <- function(n.adapt,n.chains){
  
  start.time <- Sys.time()
  cat("Model Initialized and adapted after\n")
  print(start.time)
  
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data$R <- diag(1,3,3)
  
  inits <- function(){list(x = data$y,
                           repro.mu = runif(1, 1, 3),
                           phi.l.mu = rnorm(1, 2.5, 0.1),
                           phi.n.mu = rnorm(1, 5.0631, 1.024),
                           phi.a.mu = rnorm(1, 7, 0.1),
                           grow.ln.mu = rnorm(1, -5, 0.001),
                           grow.na.mu = rbeta(1, 0.1, 1),
                           k.l2n.low = rnorm(1, 400, 10),
                           k.l2n.high = rnorm(1, 2000, 10),
                           k.n2a.low = rnorm(1, 700, 10),
                           k.n2a.high = rnorm(1, 2500, 10),
                           alpha.11 = rnorm(3,0,0.001),
                           alpha.33 = rnorm(3,0,0.001),
                           alpha.13 = rnorm(3,0,0.001),
                           alpha.21 = rnorm(3,0,0.001))}
  
  monitor <- c("x","m",
               "deviance",
               "phi.l.mu",
               "phi.n.mu",
               "phi.a.mu",
               "grow.ln.mu",
               "grow.na.mu",
               "repro.mu",
               "SIGMA",
               "theta.larvae",
               "theta.nymph",
               "theta.adult",
               "k.l2n.low", 
               "k.l2n.high",
               "k.n2a.low", 
               "k.n2a.high",
               "alpha.11",
               "alpha.33",
               "alpha.13",
               "alpha.21",
               "tau.11",
               "tau.33",
               "tau.13",
               "tau.21")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(3.7815,1.7488)     # larvae survival
  phi.n.mu ~ dnorm(3.6241,0.5005)     # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-5,1)           # larvae -> nymph transition 
  grow.na.mu ~ dbeta(0.1,1)           # nymph -> adult transiti
  repro.mu ~ dnorm(2,1)               # adult -> larvae transition (reproduction)
  k.l2n.low  ~ dnorm(400,0.0001) T(0,)
  k.l2n.high ~ dnorm(2000,0.0001) T(0,)
  k.n2a.low  ~ dnorm(700,0.0001) T(0,)
  k.n2a.high ~ dnorm(2500,0.0001) T(0,)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)              # mvn [3 x 3] site process
  tau.11 ~ dgamma(0.01,0.01)       # random site effect: larvae survival
  tau.33 ~ dgamma(0.01,0.01)       # random site effect: adult survival
  tau.13 ~ dgamma(0.01,0.01)       # random site effect: reproduction
  tau.21 ~ dgamma(0.01,0.01)       # random site effect: larvae-to-nymph
  
  ### Random site affect priors
  for(s in 1:3){
    alpha.11[s] ~ dnorm(0, tau.11)
    alpha.33[s] ~ dnorm(0, tau.33)
    alpha.13[s] ~ dnorm(0, tau.13)
    alpha.21[s] ~ dnorm(0, tau.21)
  }
  
  ### first latent process
  for(s in 1:3){
    for(l in 1:3){
      x[l, 1, s] ~ dpois(2)
    } # l
  } # s
  
  ### define parameters

  # daily nymph survival
  logit(phi.22) <- phi.n.mu 

  for(s in 1:3){
  
    # daily larvae survival by site
    logit(phi.11[s]) <- phi.l.mu + alpha.11[s]

    for(t in 1:N_days[s]){   # loop over every day in time series
      
      # larvae-to-nymph transition by site
      logit(t21[s,t]) <- grow.ln.mu + alpha.21[s]

      theta.21[s,t] <- ifelse((gdd[s,t] >= k.l2n.low) && (gdd[s,t] <= k.l2n.high),t21[s,t],0)
      theta.32[s,t] <- ifelse((gdd[s,t] <= k.n2a.low) || (gdd[s,t] >= k.n2a.high),grow.na.mu,0)
      
      A.day[1,1,s,t] <- phi.11[s]*(1-theta.21[s,t]) 
      A.day[2,1,s,t] <- phi.11[s]*theta.21[s,t] 
      A.day[2,2,s,t] <- phi.22*(1-theta.32[s,t]) 
      A.day[3,2,s,t] <- phi.22*theta.32[s,t]
      log(A.day[1,3,s,t]) <- repro.mu + alpha.13[s]
      logit(A.day[3,3,s,t]) <- phi.a.mu + alpha.33[s]
      A.day[1,2,s,t] <- 0
      A.day[2,3,s,t] <- 0
      A.day[3,1,s,t] <- 0
    }
  }
  
  ### aggregate daily matrix between sampling events
  for(s in 1:3){
    for(t in 1:(N_est[s]-1)){  ## number of days to estimate latent state
    
      ## df is the number of days between each sampling occasion
      ## dt.index the number of days elaplsed from first sampling occasion
      
      for(i in 1){  
        TRANS[1:3,1:3,s,dt.index[s,t]-df[s,t]+1] <- A.day[1:3,1:3,s,dt.index[s,t]-df[s,t]+1] 
      }                                                          
      for(i in 2:df[s,t]){
        TRANS[1:3,1:3,s,dt.index[s,t]-df[s,t]+i] <- TRANS[1:3,1:3,s,dt.index[s,t]-df[s,t]+i-1] %*% 
        A.day[1:3,1:3,s,dt.index[s,t]-df[s,t]+i]
      }
    }
  }
  
  
  ## prior probablity we observe ticks by life stage
  theta.larvae ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  ### Process Model
  for(s in 1:3){
    for(t in 1:(N_est[s]-1)){  # number of days for each site
    
    ## TRANS is the aggrigated transition matrix between observations
    
    Ex[1:3,t,s] <- TRANS[1:3,1:3,s,dt.index[s,t]] %*% x[1:3,t,s]
    x[1:3,t+1,s] ~ dmnorm(Ex[1:3,t,s], SIGMA)
    
    ## Data Model ##
    
    ## fit the blended model to observed data 
    y[1,t,s] ~ dpois(m[1,t,s])
    y[2,t,s] ~ dpois(m[2,t,s])
    y[3,t,s] ~ dpois(m[3,t,s])
    
    ## blend the poisson and zero inflation models
    m[1,t,s] <- x[1,t,s]*b.larvae[t,s] + 1E-10
    m[2,t,s] <- x[2,t,s]*b.nymph[t,s] + 1E-10
    m[3,t,s] <- x[3,t,s]*b.adult[t,s] + 1E-10
    
    ## binary outcome of observation by life stage
    b.larvae[t,s] ~ dbern(theta.larvae)
    b.nymph[t,s] ~ dbern(theta.nymph)
    b.adult[t,s] ~ dbern(theta.adult)
    
    } # t
  }
  
}" # end model

  
  load.module("glm")
  load.module("dic")
  
  j.model <- jags.model(file = textConnection(model),
                        data = data,
                        inits = inits,
                        n.adapt = n.adapt,
                        n.chains = n.chains)
  
  return <- list(j.model = j.model,
                 monitor = monitor,
                 start.time = start.time)
  
  time <- Sys.time() - start.time 
  cat("Model Initialized and adapted after\n")
  print(time)
  
  return(return)
}


