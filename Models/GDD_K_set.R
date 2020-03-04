##   This is the model for Cary Ticks on a Daily scale.    
##
##   This model has interacting parameters in the matrix (phi*(1-psi))
##   
##   Transition is a threshold (on if within gdd window, off otherwise)
##
##   No Drivers on survival or observation        
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

source("Functions/cary_tick_met_JAGS.R") # get data
source("Functions/site_data_met.R") # subset data for independent fits
source("Functions/RunMCMC_Model.R") # run mcmc and check for convergence / burnin
source("Functions/get_survival.R") # run mcmc and check for convergence / burnin
# source("Functions/RunMCMC_runjags.R") # run mcmc and check for convergence / burnin

run_model <- function(site.run, met.proc, n.adapt, n.chains){
  
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data <- site_data_met(site = site.run, met.variable = met.proc, data)
  
  seq.days <- matrix(NA, data$N_est-1, max(data$df, na.rm = TRUE))
  for(i in 1:(data$N_est-1)){
    xx <- (data$dt.index[i+1]-1):data$dt.index[i]
    seq.days[i,1:length(xx)] <- xx
  }
  data$seq.days <- seq.days
  
  data$R <- diag(1, 3, 3)
 
  # get survival estimates
  survival <- get_survival(larva.driver = met.proc,
                           nymph.driver = met.proc)
  data$larva.mean <- survival$larva.survival.mean
  data$larva.prec <- survival$larva.survival.prec
  data$nymph.mean <- survival$nymph.survival.mean
  data$nymph.prec <- survival$nymph.survival.prec
  
  inits <- function(){list(p = data$y[,-1],
                           repro.mu = runif(1, 10, 20),
                           phi.l.mu = rnorm(1, data$larva.mean, 0.1),
                           phi.n.mu = rnorm(1, data$nymph.mean, 0.1),
                           phi.a.mu = rnorm(1, 6, 0.001),
                           grow.ln.mu = rnorm(1, -6, 0.1),
                           grow.na.mu = rnorm(1, -6, 0.1))}
  
  monitor <- c("x",
               "deviance",
               "phi.l.mu",
               "phi.n.mu",
               "phi.a.mu",
               "grow.ln.mu",
               "grow.na.mu",
               "repro.mu",
               "SIGMA",
               # "theta.larva",
               "theta.nymph",
               "theta.adult",
               "beta.l.obs")
               # "beta.n.obs",
               # "beta.a.obs",
               # "beta.l.vert",
               # "beta.n.vert",
               # "beta.a.vert",
               # "beta.l.lat")
               # "beta.n.lat")
  # "beta.a.lat")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(larva.mean,larva.prec)               # larvae survival
  phi.n.mu ~ dnorm(nymph.mean,nymph.prec)               # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-6, 1)             # larvae -> nymph transition 
  grow.na.mu ~ dnorm(-6, 1)             # nymph -> adult transition
  repro.mu ~ dnorm(7, 0.01) T(0,)              # adult -> larvae transition (reproduction)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process
  
  ## observation regression priors
  beta.l.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.n.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.a.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.l.vert ~ dnorm(0, 0.001) T(0,)
  # beta.n.vert ~ dnorm(0, 0.001) T(0,)
  # beta.a.vert ~ dnorm(0, 0.001) T(0,)
  # beta.l.lat ~ dnorm(0, 0.001)
  # beta.n.lat ~ dnorm(0, 0.001)
  # beta.a.lat ~ dnorm(0, 0.001)
  # theta.larva ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  ### first latent process 
  x[1, 1] ~ dpois(1) 
  x[2, 1] ~ dpois(10) 
  x[3, 1] ~ dpois(1) 
  
  ## missing temperature model - observation
  for(t in met.obs.miss){
  met.obs[t] ~ dunif(met.obs.range[1], met.obs.range[2])
  }
  
  logit(phi.11) <- phi.l.mu 
  logit(phi.22) <- phi.n.mu
  logit(l2n) <- grow.ln.mu
  logit(n2a) <- grow.na.mu
  
  ### define parameters
  for(t in 1:N_days){   # loop over every day in time series
  
  ## Survival parameters are random intercept plus fixed effect on temperature
  ## transition is a threshold (on if within gdd window, off otherwise)
  
  theta.21[t] <- ifelse((gdd[t] >= 500) && (gdd[t] <= 2500),l2n,0)
  theta.32[t] <- ifelse((gdd[t] <= 1000) || (gdd[t] >= 2500),n2a,0)
  lambda[t] <- ifelse((gdd[t] >= 1500) && (gdd[t] <= 2500),repro.mu,0)
  
  A.day[1,1,t] <- phi.11*(1-theta.21[t]) 
  A.day[2,1,t] <- phi.11*theta.21[t] 
  A.day[2,2,t] <- phi.22*(1-theta.32[t]) 
  A.day[3,2,t] <- phi.22*theta.32[t]
  logit(A.day[3,3,t]) <- phi.a.mu 
  A.day[1,3,t] <- lambda[t]
  A.day[1,2,t] <- 0
  A.day[2,3,t] <- 0
  A.day[3,1,t] <- 0
  }
  
  ### aggregate daily matrix between sampling events
  
  for(i in 1:(N_est-1)){  ## number of days to estimate latent state
  
    for(day in seq.days[i,1]){
      TRANS[1:3,1:3,day] <- A.day[,,day] %*% A.day[,,day-1]
    } 

    for(day in seq.days[i,2:df[i]]){
      TRANS[1:3,1:3,day] <- TRANS[1:3,1:3,day+1] %*% A.day[,,day]
    }
  }
  
  
  ### Process Model
  
  for(t in 1:(N_est-1)){
  
  # expected number questing
  Ex[1:3,t] <- TRANS[1:3,1:3,dt.index[t]] %*% x[1:3,t]
  
  # process error
  p[1:3,t] ~ dmnorm(Ex[1:3,t], SIGMA)
  x[1,t+1] <- max(p[1,t], 0)
  x[2,t+1] <- max(p[2,t], 0)
  x[3,t+1] <- max(p[3,t], 0)
  }
  
  ### Data Model ###
  for(t in 1:N_est){
  
  ## fit the blended model to observed data 
  y[1,t] ~ dpois(m[1,t])
  y[2,t] ~ dpois(m[2,t])
  y[3,t] ~ dpois(m[3,t])
  
  ## blend the poisson and zero inflation models
  m[1,t] <- x[1,t]*b.larva[t] + 1E-10
  m[2,t] <- x[2,t]*b.nymph[t] + 1E-10
  m[3,t] <- x[3,t]*b.adult[t] + 1E-10
  
  ## observation probability based on temperature
  theta.larva[t] <- 1 / (1 + beta.l.obs*(met.obs[t])^2)
  # theta.nymph[t] <- 1 / (1 + beta.n.obs*(met.obs[t] + beta.n.lat)^2)
  # theta.adult[t] <- 1 / (1 + beta.a.vert + beta.a.obs*(met.obs[t] + beta.a.lat)^2)
  
  ## binary outcome of observation by life stage
  b.larva[t] ~ dbern(theta.larva[t])
  b.nymph[t] ~ dbern(theta.nymph)
  b.adult[t] ~ dbern(theta.adult)
  
  } # t
}" 
  
  load.module("glm")
  load.module("dic")
  
  start.time <- Sys.time()
  
  cat("\nCompiling model with", n.adapt, "adaptation iterations\n")
  
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



# 
# print("-------------- DIC ----------------")
# 
# dic.samples(j.model, 10000)




