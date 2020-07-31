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

run_model <- function(site.run, met.proc, n.adapt, n.chains) {
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data <- site_data_met(site = site.run, 
                        met.variable = met.proc, 
                        data = data,
                        time.effect = "month")
  
  # met.diff <- rep(NA, data$N_days)
  # for(m in 2:length(data$met)){
  #   met.diff[m] <- data$met[m] - data$met[m-1]
  # }
  # data$met.diff <- met.diff
  # data$met.mis.diff <- which(is.na(met.diff))
  # data$met.range.diff <- range(met.diff, na.rm = TRUE)
  
  # which transitions are set and which are variable
  data$rho.n.mu <- 500
  data$rho.l.mu <- 1500
  data$rho.a.mu <- 2500
  data$rho.prec <- 1 / 200^2
  
  # get survival estimates
  survival <- get_survival(larva.driver = met.proc,
                           nymph.driver = NULL)
  data$larva.mean <- survival$larva.survival.mean
  data$larva.prec <- survival$larva.survival.prec
  data$nymph.mean <- survival$nymph.survival.mean
  data$nymph.prec <- survival$nymph.survival.prec
  
  data$larva.beta.mu <- survival$larva.beta.mean
  data$larva.beta.prec <- survival$larva.beta.prec
  
  inits <- function() {
    list(
      p = data$y[, -1],
      repro.mu = runif(1, 0, 10),
      phi.l.mu = rnorm(1, data$larva.mean, 0.1),
      phi.n.mu = rnorm(1, data$nymph.mean, 0.1),
      phi.a.mu = rnorm(1, 6, 0.001),
      rho.a = runif(1, 2400, 2600),
      rho.l = runif(1, 1400, 1600),
      rho.n = runif(1, 400, 600),
      grow.ln.mu = rnorm(1, -6, 0.1),
      grow.na.mu = rnorm(1, -6, 0.1)
    )
  }
  
  monitor <- c(
    "x",
    "deviance",
    "phi.l.mu",
    "phi.n.mu",
    "phi.a.mu",
    "beta.l",
    "grow.ln.mu",
    "grow.na.mu",
    "repro.mu",
    "SIGMA",
    "theta.nymph",
    "theta.adult",
    "beta.l.obs",
    "rho.l",
    "rho.n",
    "rho.a",
    "alpha.month"
  )

  model = " model {

  ### global priors
  phi.l.mu ~ dnorm(larva.mean,larva.prec)               # larvae survival
  phi.n.mu ~ dnorm(nymph.mean,nymph.prec)               # nymph survival
  phi.a.mu ~ dnorm(5,1)               # adult survival
  grow.ln.mu ~ dnorm(-6, 1)             # larvae -> nymph transition
  grow.na.mu ~ dnorm(-6, 1)             # nymph -> adult transition
  repro.mu ~ dnorm(7, 0.01) T(0,)              # adult -> larvae transition (reproduction)
  rho.l ~ dnorm(rho.l.mu, rho.prec)
  rho.n ~ dnorm(rho.n.mu, rho.prec)
  rho.a ~ dnorm(rho.a.mu, rho.prec)

  ### precision priors
  SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process
  
  ### random month prior
  for(s in 1:3){
    for(month in n.months){
      alpha.month[s,month] ~ dnorm(0, 0.001)
    }
  }
  
  beta.l ~ dnorm(larva.beta.mu, larva.beta.prec)

  ## observation regression priors
  beta.l.obs ~ dnorm(0, 0.001) T(1E-10,)
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

  ## missing process met
  for(t in met.mis){
    met[t] ~ dunif(met.range[1], met.range[2])
  }
  
  # for(t in met.mis.diff){
  #   met.diff[t] ~ dunif(met.range.diff[1], met.range.diff[2])
  # }
  
  logit(phi.22) <- phi.n.mu
  logit(l2n) <- grow.ln.mu
  logit(n2a) <- grow.na.mu

  ### define parameters
  for(t in 1:N_days){   # loop over every day in time series

    theta.21[t] <- ifelse((gdd[t] >= rho.n) && (gdd[t] <= 2500),l2n,0)
    theta.32[t] <- ifelse((gdd[t] <= 1000) || (gdd[t] >= rho.a),n2a,0)
    lambda[t] <- ifelse((gdd[t] >= rho.l) && (gdd[t] <= 2500),repro.mu,0)
    logit(phi.11[t]) <- phi.l.mu + beta.l*met[t]
  
    A.day[1,1,t] <- phi.11[t]*(1-theta.21[t])
    A.day[2,1,t] <- phi.11[t]*theta.21[t]
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
    Ex[1:3,t] <- TRANS[1:3,1:3,dt.index[t]] %*% x[1:3,t] + alpha.month[1:3, month.index[t]]
  
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
  
  j.model <- jags.model(
    file = textConnection(model),
    data = data,
    inits = inits,
    n.adapt = n.adapt,
    n.chains = n.chains
  )
  
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
