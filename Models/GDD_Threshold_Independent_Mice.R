##   This is the model for Cary Ticks on a Daily scale.    
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

#'@param site.run Which site is running? One of "Green Control", "Henry Conrtol", "Tea Control"
#'@param n.adapt number of adaption iterations
#'@param n.chains number of chains, 1 for parallel batch jobs
#'@param burnin number of burnin iterations
#'@param thin thinning interval
#'@param n.iter number of iterations post burnin.     

library(ecoforecastR)

source("Functions/site_data_met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/mice_estimated_jags.R")

run_model <- function(site.run,n.adapt,n.chains){
  
  start.time <- Sys.time()
  cat("Model Initialized and adapted after\n")
  print(start.time)
  
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data <- site_data_met(site = site.run, met.variable = NULL, data)
  data$R <- diag(1,3,3)
  
  # get mice data
  mice <- mice_estimated_jags(site.run)
  data$mice.mean <- mice$mice.mean 
  data$mice.prec <- mice$mice.prec 
    
    
inits <- function(){list(x = data$y,
                         repro.mu = runif(1, 1, 3),
                         phi.l.mu = rnorm(1, 3.501, 1.4993),
                         phi.n.mu = rnorm(1, 5.0631, 1.024),
                         phi.a.mu = rnorm(1, 3, 1),
                         grow.ln.mu = rnorm(1, -5, 0.001),
                         grow.na.mu = rnorm(1, -5, 0.001),
                         beta.21 = rnorm(1, 0, 0.001),
                         beta.32 = rnorm(1, 0, 0.001))}

monitor <- c("x","m",
             "deviance",
             "phi.l.mu",
             "phi.n.mu",
             "phi.a.mu",
             "grow.ln.mu",
             "grow.na.mu",
             "repro.mu",
             "SIGMA",
             "beta.21",
             "beta.32",
             "mice",
             "theta.larvae",
             "theta.nymph",
             "theta.adult")

model = " model {

### global hyperpriors
phi.l.mu ~ dnorm(3.7815,1.7488)               # larvae survival
phi.n.mu ~ dnorm(3.6241,0.5005)               # nymph survival
phi.a.mu ~ dnorm(5,1)               # adult survival
grow.ln.mu ~ dnorm(-5,0.01)             # larvae -> nymph transition 
grow.na.mu ~ dnorm(-5,0.01)             # nymph -> adult transiti
repro.mu ~ dnorm(2,1)              # adult -> larvae transition (reproduction)

### precision priors
SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process

### beta priors
beta.21 ~ dnorm(0, 0.1)
beta.32 ~ dnorm(0, 0.1)

### first latent process 
for(l in 1:3){
  x[l, 1] ~ dpois(2) 
} # l

logit(phi.11) <- phi.l.mu 
logit(phi.22) <- phi.n.mu 

### define parameters
for(t in 1:N_days){   # loop over every day in time series
  
  ### EIV mice
  mice[t] ~ dnorm(mice.mean[t], mice.prec[t]) T(0,)

  logit(t21[t]) <- grow.ln.mu + beta.21*mice[t]
  logit(t32[t]) <- grow.na.mu + beta.32*mice[t]

  # theta.21[t] <- ifelse((gdd[t] >= 500),t21[t],0)
  theta.21[t] <- ifelse((gdd[t] >= 500) && (gdd[t] <= 2000),t21[t],0)

  theta.32[t] <- ifelse((gdd[t] <= 750) || (gdd[t] >= 2500),t32[t],0)

  A.day[1,1,t] <- phi.11*(1-theta.21[t]) 
  A.day[2,1,t] <- phi.11*theta.21[t] 
  A.day[2,2,t] <- phi.22*(1-theta.32[t]) 
  A.day[3,2,t] <- phi.22*theta.32[t]
  logit(A.day[3,3,t]) <- phi.a.mu 
  log(A.day[1,3,t]) <- repro.mu
  A.day[1,2,t] <- 0
  A.day[2,3,t] <- 0
  A.day[3,1,t] <- 0
}

### aggregate daily matrix between sampling events

for(t in 1:(N_est-1)){  ## number of days to estimate latent state

  ## df is the number of days between each sampling occasion
  ## dt.index the number of days elaplsed from first sampling occasion
  
  for(i in 1){  
    TRANS[1:3,1:3,dt.index[t]-df[t]+1] <- A.day[1:3,1:3,dt.index[t]-df[t]+1] 
  }                                                          
  for(i in 2:df[t]){
    TRANS[1:3,1:3,dt.index[t]-df[t]+i] <- TRANS[1:3,1:3,dt.index[t]-df[t]+i-1] %*% 
    A.day[1:3,1:3,dt.index[t]-df[t]+i]
  }
}


## prior probablity we observe ticks by life stage
theta.larvae ~ dunif(0,1)
theta.nymph ~ dunif(0,1)
theta.adult ~ dunif(0,1)

### Process Model
for(t in 1:(N_est-1)){  # number of days for each site
  
    ## TRANS is the aggrigated transition matrix between observations
    
  Ex[1:3,t] <- TRANS[1:3,1:3,dt.index[t]] %*% x[1:3,t]
  x[1:3,t+1] ~ dmnorm(Ex[1:3,t], SIGMA)
    
    ## Data Model ##

    ## fit the blended model to observed data 
    y[1,t] ~ dpois(m[1,t])
    y[2,t] ~ dpois(m[2,t])
    y[3,t] ~ dpois(m[3,t])
    
    ## blend the poisson and zero inflation models
    m[1,t] <- x[1,t]*b.larvae[t] + 1E-10
    m[2,t] <- x[2,t]*b.nymph[t] + 1E-10
    m[3,t] <- x[3,t]*b.adult[t] + 1E-10
    
    ## binary outcome of observation by life stage
    b.larvae[t] ~ dbern(theta.larvae)
    b.nymph[t] ~ dbern(theta.nymph)
    b.adult[t] ~ dbern(theta.adult)

  } # t

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


