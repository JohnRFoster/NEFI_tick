##   This is the model for Cary Ticks on a Daily scale.    
##
##   This model has interactinf parameters in the matrix (phi*(1-psi))
##   
##   Transition is a threshold (on if within gdd window, off otherwise)
##
##   Survival in the questing matrix is estimated with the a random intercept, and a 
##   fixed effect on temperature.        
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

run_model <- function(site.run, n.adapt, n.chains){
  
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data <- site_data_met(site = site.run, met.variable = "temp", data)
  data$R <- diag(1,3,3)
  
  inits <- function(){list(p = data$y[,-1],
                           repro.mu = rnorm(1, 3, 0.001),
                           phi.l.mu = rnorm(1, 2, 0.001),
                           phi.n.mu = rnorm(1, 4, 0.001),
                           phi.a.mu = rnorm(1, 6, 0.001),
                           grow.ln.mu = rbeta(1, 0.1, 1),
                           grow.na.mu = rbeta(1, 0.1, 1),
                           beta.11 = rnorm(1, 0.11420, 0.001),
                           beta.22 = rnorm(1, -0.05054, 0.001),
                           beta.33 = runif(1, 0, 2))}
  
  monitor <- c("x","m",
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
               "k.21.low",
               "k.21.high",
               "k.32.low",
               "k.32.high",
               "k.13.low",
               "k.13.high",
               "theta.larvae",
               "theta.nymph",
               "theta.adult")
  
  model = " model {
  
  ### global hyperpriors
  phi.l.mu ~ dnorm(3.07431,1.2228)               # larvae survival
  phi.n.mu ~ dnorm(4.15213,1.5276)               # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dbeta(0.1,1)             # larvae -> nymph transition 
  grow.na.mu ~ dbeta(0.1,1)             # nymph -> adult transition
  repro.mu ~ dnorm(2, 1)              # adult -> larvae transition (reproduction)
  
  ### precision priors
  SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process
  
  ### Fixed effects priors
  beta.11 ~ dnorm(0.1142,2104.1999)
  beta.22 ~ dnorm(-0.05054,383.5644)
  beta.33 ~ dgamma(0.9,1)
  k.21.low  ~  dnorm(500,0.001) T(0,)
  k.21.high ~ dnorm(2500,0.001) T(0,)
  k.32.low  ~  dnorm(1000,0.001) T(0,)
  k.32.high ~ dnorm(2500,0.001) T(0,)
  k.13.low  ~  dnorm(1500,0.001) T(0,)
  k.13.high ~ dnorm(2500,0.001) T(0,)

  ## prior probablity we observe ticks by life stage
  theta.larvae ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  for(t in 1:39){
    met[temp.mis[t]] ~ dunif(-30,23)
  }
  
  ### first latent process 
  for(l in 1:3){
    x[l, 1] ~ dpois(10) 
  } # l
  
  ### define parameters
  for(t in 1:N_days){   # loop over every day in time series
  
    ## Survival parameters are random intercept plus fixed effect on temperature
    ## transition is a threshold (on if within gdd window, off otherwise)

    logit(phi.11[t]) <- phi.l.mu + beta.11*met[t]
    logit(phi.22[t]) <- phi.n.mu + beta.22*met[t]
    
    theta.21[t] <- ifelse((gdd[t] >= k.21.low) && (gdd[t] <= k.21.high),grow.ln.mu,0)
    theta.32[t] <- ifelse((gdd[t] <= k.32.high) || (gdd[t] >= k.32.low),grow.na.mu,0)
    log(lambda[t]) <- ifelse((gdd[t] >= k.13.low) && (gdd[t] <= k.13.high),repro.mu,0)

    A.day[1,1,t] <- phi.11[t]*(1-theta.21[t]) 
    A.day[2,1,t] <- phi.11[t]*theta.21[t] 
    A.day[2,2,t] <- phi.22[t]*(1-theta.32[t]) 
    A.day[3,2,t] <- phi.22[t]*theta.32[t]
    logit(A.day[3,3,t]) <- phi.a.mu + beta.33*met[t]
    A.day[1,3,t] <- lambda[t]
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
  
  ### Process Model
  
  for(t in 1:(N_est-1)){

  # expected number questing
  Ex[1:3,t] <- TRANS[1:3,1:3,dt.index[t]] %*% x[1:3,t]

  # process error
  p[1:3,t] ~ dmnorm(Ex[1:3,t], SIGMA)
  x[1,t+1] <- max(p[1,t], 0)
  x[2,t+1] <- max(p[2,t], 0)
  x[3,t+1] <- max(p[3,t], 0)
  
  ### Data Model ###
  
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
}" 
  
  load.module("glm")
  load.module("dic")
  
  start.time <- Sys.time()
  
  cat("\nCompiling model with", n.adapt, "adaptation iterations \n")

  j.model <- jags.model(file = textConnection(model),
                        data = data,
                        inits = inits,
                        n.adapt = n.adapt,
                        n.chains = n.chains)
  
  time <- Sys.time() - start.time 
  cat("Model Initialized and adapted after\n")
  print(time)
  
  jags.out <- runMCMC_Model(j.model, variableNames = monitor)
  
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



# 
# print("-------------- DIC ----------------")
# 
# dic.samples(j.model, 10000)




