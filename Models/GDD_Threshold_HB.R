##   This is the HB model for Cary Ticks on a Daily scale.  
##
##   Random site effects on nymph survival, larvae-to-nymph transition, and reproduction
##
##   This model has interacting parameters in the matrix (phi*(1-psi))
##   
##   Transition is a threshold (on if within gdd window, off otherwise)
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
run_model <- function(n.adapt,n.chains,burnin,thin,n.iter){
  
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
                           alpha.22 = rnorm(3,0,0.1),
                           alpha.13 = rnorm(3,0,0.1),
                           alpha.32 = rnorm(3,0,0.1))}
  
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
               "alpha.22",
               "alpha.21",
               "alpha.13",
               "tau.22",
               "tau.21",
               "tau.13")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(3.7815,1.7488)     # larvae survival
  phi.n.mu ~ dnorm(3.6241,0.5005)     # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-5,1)           # larvae -> nymph transition 
  grow.na.mu ~ dbeta(0.1,1)           # nymph -> adult transiti
  repro.mu ~ dnorm(2,1)               # adult -> larvae transition (reproduction)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)              # mvn [3 x 3] site process
  tau.22 ~ dgamma(0.01,0.01)       # random site effect: nymph survival
  tau.21 ~ dgamma(0.01,0.01)       # random site effect: larvae-adult transition
  tau.13 ~ dgamma(0.01,0.01)       # random site effect: reproduction

  ### Random site affect priors
  for(s in 1:3){
    alpha.22[s] ~ dnorm(0, tau.22)
    alpha.21[s] ~ dnorm(0, tau.21)
    alpha.13[s] ~ dnorm(0, tau.13)
  }
  
  ### first latent process
  for(s in 1:3){
    for(l in 1:3){
      x[l, 1, s] ~ dpois(2)
    } # l
  } # s
  
  ### define parameters
  for(s in 1:3){
    
    ## daily nymph survival by site
    logit(phi.22[s]) <- phi.n.mu + alpha.22[s]
    
    for(t in 1:N_days[s]){   # loop over every day in time series
    
      ## Survival parameters are random intercept plus fixed effect on temperature
      ## transition is a threshold (on if within gdd window, off otherwise)

      logit(t21[s,t]) <- grow.ln.mu + alpha.21[s]
      theta.21[s,t] <- ifelse((gdd[s,t] >= 500),t21[s,t],0)
      theta.32[s,t] <- ifelse((gdd[s,t] <= 750) || (gdd[s,t] >= 2500),grow.na.mu,0)
      
      A.day[1,1,s,t] <- phi.11*(1-theta.21[s,t]) 
      A.day[2,1,s,t] <- phi.11*theta.21[s,t] 
      A.day[2,2,s,t] <- phi.22[s]*(1-theta.32[s,t]) 
      A.day[3,2,s,t] <- phi.22[s]*theta.32[s,t]
      log(A.day[1,3,s,t]) <- repro.mu + alpha.13[s]
      logit(A.day[3,3,s,t]) <- phi.a.mu
      A.day[1,2,s,t] <- 0
      A.day[2,3,s,t] <- 0
      A.day[3,1,s,t] <- 0
    }
  }
  
  logit(phi.11) <- phi.l.mu 
  
  
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
  
  time <- Sys.time() - start.time 
  cat("Model Initialized and adapted after\n")
  print(time)
  
  jags.out <- coda.samples(model = j.model,
                           variable.names = monitor,
                           burnin = burnin,
                           thin = thin,
                           n.iter = n.iter)
  
  time <- Sys.time() - start.time 
  cat(n.iter,"iterations complete after\n")
  print(time)
  
  ## split output
  out <- list(params = NULL, predict = NULL, m.cols = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  m.cols <- grep("m[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$m.cols <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, m.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -c(m.cols, pred.cols)])
  out <- list(out = out, j.model = j.model)
  
  cat("\nJAGS model split and returned\n")
  
  time <- Sys.time() - start.time 
  cat("\nTotal time\n")
  print(time)
  
  return(out)
  }


