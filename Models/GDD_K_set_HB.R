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
source("Functions/hb_restructure.R")
source("Functions/get_inits_HB.R")

run_model <- function(met.proc, n.adapt, n.chains){
  
  data <- cary_ticks_met_JAGS()
  data <- hb_restructure(data)
  data$R <- diag(1, 3, 3)
  
  # get survival estimates
  survival <- get_survival(larva.driver = met.proc,
                           nymph.driver = met.proc)
  data$larva.mean <- survival$larva.survival.mean
  data$larva.prec <- survival$larva.survival.prec
  data$nymph.mean <- survival$nymph.survival.mean
  data$nymph.prec <- survival$nymph.survival.prec
  
  # get inits from independent models
  init.params <- get_inits_HB(met.proc)
  
  inits <- function(){list(p = data$y[,-1,],
                           repro.mu = rnorm(1, init.params["repro.mu"], 1),
                           phi.l.mu = rnorm(1, init.params["phi.l.mu"], 0.01),
                           phi.n.mu = rnorm(1, init.params["phi.n.mu"], 0.01),
                           phi.a.mu = rnorm(1, init.params["phi.a.mu"], 0.01),
                           grow.ln.mu = rnorm(1, init.params["grow.ln.mu"], 0.1),
                           grow.na.mu = rnorm(1, init.params["grow.na.mu"], 0.1),
                           theta.larva = runif(1, init.params["theta.larva"]*0.9, 0.9999),
                           theta.nymph = runif(1, init.params["theta.nymph"]*0.9, 0.9999),
                           theta.adult = runif(1, init.params["theta.adult"]*0.9, 0.9999),
                           alpha.a = rnorm(3, 0, 1))}
  
  data$temp.min <- data$temp.max <- data$rh.max <- data$rh.min <- data$precip <- data$vpd <- NULL
  data$mis.rh.max <- data$mis.rh.min <- data$mis.temp.max <- data$mis.temp.min <- data$mis.vpd <- NULL
  
  monitor <- c("x",
               "deviance",
               "phi.l.mu",
               "phi.n.mu",
               "phi.a.mu",
               "grow.ln.mu",
               "grow.na.mu",
               "repro.mu",
               "SIGMA",
               "theta.larva",
               "theta.nymph",
               "theta.adult",
               # "beta.l.obs",
               # "beta.l.lat",
               # "beta.l.vert")
               # "beta.n.obs",
               # "beta.n.lat",
               # "beta.n.vert",
               # "beta.a.obs",
               # "beta.a.lat",
               # "beta.a.vert")
               "alpha.a")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(larva.mean,larva.prec)               # larvae survival
  phi.n.mu ~ dnorm(nymph.mean,nymph.prec)               # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-6, 1)             # larvae -> nymph transition 
  grow.na.mu ~ dnorm(-6, 1)             # nymph -> adult transition
  repro.mu ~ dnorm(7, 0.01) T(0,)              # adult -> larvae transition (reproduction)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)          # mvn [3 x 3] site process
  tau.a ~ dgamma(0.001, 0.001) # adult survival random effect 

  ### random effect prior
  for(s in 1:3){
    alpha.a[s] ~ dnorm(0, tau.a)
  }
  
  ## observation regression priors
  # beta.l.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.n.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.a.obs ~ dnorm(0, 0.001) T(1E-10,)
  # beta.l.vert ~ dnorm(0, 0.001) T(0,)
  # beta.n.vert ~ dnorm(0, 0.001) T(0,)
  # beta.a.vert ~ dnorm(0, 0.001) T(0,)
  # beta.l.lat ~ dnorm(0, 0.001)
  # beta.n.lat ~ dnorm(0, 0.001)
  # beta.a.lat ~ dnorm(0, 0.001)
  theta.larva ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  ### first latent process 
  for(s in 1:3){
    x[1, 1, s] ~ dpois(1) 
    x[2, 1, s] ~ dpois(10) 
    x[3, 1, s] ~ dpois(1) 
  }
  
  logit(phi.11) <- phi.l.mu 
  logit(phi.22) <- phi.n.mu
  logit(l2n) <- grow.ln.mu
  logit(n2a) <- grow.na.mu
  
  
    ### define parameters
  for(t in 1:max(N_days)){   # loop over every day in time series
  
    ## Survival parameters are random intercept plus fixed effect on temperature
    ## transition is a threshold (on if within gdd window, off otherwise)
    
    theta.21[t] <- ifelse((gdd[t] >= 500) && (gdd[t] <= 2500),l2n,0)
    theta.32[t] <- ifelse((gdd[t] <= 1000) || (gdd[t] >= 2500),n2a,0)
    lambda[t] <- ifelse((gdd[t] >= 1500) && (gdd[t] <= 2500),repro.mu,0)
      
    for(s in 1:3){
      A.day[1,1,t,s] <- phi.11*(1-theta.21[t]) 
      A.day[2,1,t,s] <- phi.11*theta.21[t] 
      A.day[2,2,t,s] <- phi.22*(1-theta.32[t]) 
      A.day[3,2,t,s] <- phi.22*theta.32[t]
      logit(A.day[3,3,t,s]) <- phi.a.mu + alpha.a[s]
      A.day[1,3,t,s] <- lambda[t]
      A.day[1,2,t,s] <- 0
      A.day[2,3,t,s] <- 0
      A.day[3,1,t,s] <- 0
    }
  }
  
  
  ### aggregate daily matrix between sampling events
  for(s in 1:3){
    for(i in 1:(N_est[s]-1)){  ## number of days to estimate latent state
    
      for(day in seq.days[i,1,s]){
        TRANS[1:3,1:3,day,s] <- A.day[,,day,s] %*% A.day[,,day-1,s]
      } 
  
      for(day in seq.days[i,2:df[s,i],s]){
        TRANS[1:3,1:3,day,s] <- TRANS[1:3,1:3,day+1,s] %*% A.day[,,day,s]
      }
    }
  }
  
  ### Process Model
  for(s in 1:3){
    for(t in 1:(N_est[s]-1)){
  
      # expected number questing
      Ex[1:3,t,s] <- TRANS[1:3,1:3,dt.index[s,t],s] %*% x[1:3,t,s]
      
      # process error
      p[1:3,t,s] ~ dmnorm(Ex[1:3,t,s], SIGMA)
      x[1,t+1,s] <- max(p[1,t,s], 0)
      x[2,t+1,s] <- max(p[2,t,s], 0)
      x[3,t+1,s] <- max(p[3,t,s], 0)
    }
  }
  
  
  ### Data Model ###
  for(s in 1:3){
    for(t in 1:N_est[s]){
  
      ## fit the blended model to observed data 
      y[1,t,s] ~ dpois(m[1,t,s])
      y[2,t,s] ~ dpois(m[2,t,s])
      y[3,t,s] ~ dpois(m[3,t,s])
      
      ## blend the poisson and zero inflation models
      m[1,t,s] <- x[1,t,s]*b.larva[t,s] + 1E-10
      m[2,t,s] <- x[2,t,s]*b.nymph[t,s] + 1E-10
      m[3,t,s] <- x[3,t,s]*b.adult[t,s] + 1E-10
      
      ## observation probability based on temperature
      # theta.larva[t] <- 1 / (1 + beta.l.obs*(met.obs[t])^2)
      # theta.nymph[t] <- 1 / (1 + beta.n.obs*(met.obs[t])^2)
      # theta.adult[t] <- 1 / (1 + beta.a.vert + beta.a.obs*(met.obs[t])^2)
      
      ## binary outcome of observation by life stage
      b.larva[t,s] ~ dbern(theta.larva)
      b.nymph[t,s] ~ dbern(theta.nymph)
      b.adult[t,s] ~ dbern(theta.adult)
      
    } # t
  }
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




