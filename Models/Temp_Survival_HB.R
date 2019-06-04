##   This is the model for Cary Ticks on a Daily scale.    
##
##   Runs for all sites - HB
##   
##   Transition is estimated with a random intercept.
##
##   Survival is estimated with the a fixed random intercept, a random effect by site, and a
##   fixed effect on temperature
##
##   The state is estimated for every sampling day only, 
##   but demographic params are estimated daily.                       
##
##   Zero Inflated Poisson Data model by life stage.  

#'@param site.run Which site is running? One of "Green Control", "Henry Conrtol", "Tea Control"
#'@param n.adapt number of adaption iterations
#'@param n.chains number of chains, 1 for parallel batch jobs
#'@param burnin number of burnin iterations
#'@param thin thinning interval
#'@param n.iter number of iterations post burnin

run_model <- function(site.run,n.adapt,n.chains,burnin,thin,n.iter){
  
  start.time <- Sys.time()
  cat("Model Initialized and adapted after\n")
  print(start.time)
  
  sites <- c("Green Control","Henry Control","Tea Control")
  data <- cary_ticks_met_JAGS(sites)
  data$R <- diag(1,3,3)
  
  inits <- function(){list(x = data$y,
                           repro.mu = rnorm(1, 1, 0.001),
                           phi.l.mu = rnorm(1, 3.07431, 0.001),
                           phi.n.mu = rnorm(1, 4.15213, 0.001),
                           phi.a.mu = rnorm(1, 5, 0.001),
                           grow.ln.mu = rnorm(1, -3.07431, 0.001),
                           grow.na.mu = rnorm(1, -4.15213, 0.001),
                           alpha.11 = rnorm(3, 0, 0.001),
                           alpha.22 = rnorm(3, 0, 0.001),
                           alpha.33 = rnorm(3, 0, 0.001),
                           alpha.13 = rnorm(3, 0, 0.001),
                           alpha.32 = rnorm(3, 0, 0.001),
                           beta.11 = rnorm(1, 0.11420, 0.001),
                           beta.22 = rnorm(1, -0.05054, 0.001),
                           beta.33 = runif(1, 0, 10))}
  
  monitor <- c("x", 
               "deviance",
               "phi.l.mu",
               "phi.n.mu",
               "phi.a.mu",
               "grow.ln.mu",
               "grow.na.mu",
               "repro.mu",
               "SIGMA",
               "beta.11",
               "beta.22",
               "beta.33",
               "alpha.11",
               "alpha.22",
               "alpha.33",
               "alpha.13",
               "alpha.32",
               "tau.11",
               "tau.22",
               "tau.33",
               "tau.13",
               "tau.32",
               "theta.larvae",
               "theta.nymph",
               "theta.adult")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(3.07431,1.2228)               # larvae survival
  phi.n.mu ~ dnorm(4.15213,1.5276)               # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-3.07431,0.01)             # larvae -> nymph transition 
  grow.na.mu ~ dnorm(-4.15213,0.01)             # nymph -> adult transition
  repro.mu ~ dnorm(2, 1)              # adult -> larvae transition (reproduction)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process
  tau.11 ~ dgamma(0.01,0.01)       # random site effect: larva survival
  tau.22 ~ dgamma(0.01,0.01)       # random site effect: nymph survival
  tau.33 ~ dgamma(0.01,0.01)       # random site effect: adult survival
  tau.32 ~ dgamma(0.01,0.01)       # random site effect: nymph-adult survival
  tau.13 ~ dgamma(0.01,0.01)       # random site effect: reproduction
  
  ### Random site affect priors
  for(s in 1:3){
  alpha.11[s] ~ dnorm(0, tau.11)
  alpha.22[s] ~ dnorm(0, tau.22)
  alpha.33[s] ~ dnorm(0, tau.33)
  alpha.32[s] ~ dnorm(0, tau.32)
  alpha.13[s] ~ dnorm(0, tau.13)
  }
  
  ### Fixed effect - Daily RH - priors
  beta.11 ~ dnorm(0.1142,2104.1999)
  beta.22 ~ dnorm(-0.05054,383.5644)
  beta.33 ~ dnorm(0,0.01)
  
  ## prior probablity we observe ticks by life stage
  theta.larvae ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  ### Missing temp priors
  for(s in 1:3){
  for(t in 1:38){
  met[temp.mis[t,s],1,s] ~ dunif(-30,23)
  }
  }
  
  ### first latent process 
  for(s in 1:3){
  for(l in 1:3){
  x[l, 1, s] ~ dpois(2) 
  } # l
  } # s
  
  ### define parameters
  for(s in 1:3){
  for(t in 1:N_days[s]){   # loop over every day in time series
  
  ## site specific  daily transtion matrices;
  ## parameters are a linear combination of
  ## the global parameter (intercept) and a 
  ## site random effect
  
  logit(A.day[1,1,s,t]) <- phi.l.mu + beta.11*met[t,1,s] + alpha.11[s]
  logit(A.day[2,1,s,t]) <- grow.ln.mu 
  logit(A.day[2,2,s,t]) <- phi.n.mu + beta.22*met[t,1,s] + alpha.22[s]
  logit(A.day[3,2,s,t]) <- grow.na.mu + alpha.32[s] 
  logit(A.day[3,3,s,t]) <- phi.a.mu + beta.33*met[t,1,s] + alpha.33[s] 
  log(A.day[1,3,s,t]) <- repro.mu + alpha.13[s] 
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
  
  
  ### Process Model
  for(s in 1:3){
  for(t in 1:(N_est[s]-1)){
  
  ## within site process model, latent state
  ## normally distributed around expected value
  ## with whithin-site variability
  
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
  } # s
  
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
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  out <- list(out = out, j.model = j.model)
  
  cat("\nJAGS model split and returned\n")
  
  time <- Sys.time() - start.time 
  cat("\nTotal time\n")
  print(time)
  
  return(out)
}
