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
source("Functions/get_survival.R") # run mcmc and check for convergence / burnin
# source("Functions/RunMCMC_runjags.R") # run mcmc and check for convergence / burnin

run_model <- function(site.run, met.proc, n.adapt, n.chains){
  
  data <- cary_ticks_met_JAGS()
  
  # subset data to site.run and met variable of interest
  data <- site_data_met(site = site.run, met.variable = met.proc, data)
  data$R <- diag(1, 3, 3)
  data <- keep(data, names(data) %in% c("y", "N_est", "R"))
  
  inits <- function(){list(p = data$y[,-1])}
  
  monitor <- c("x",
               "deviance",
               "SIGMA")

  model = " model {
  
  # precision priors
  SIGMA ~ dwish(R, 4)         # mvn [3 x 3] site process
  
  # first latent state priors
  x[1,1] ~ dpois(1) 
  x[2,1] ~ dpois(10) 
  x[3,1] ~ dpois(1)  
  
  # data model
  for(t in 1:N_est){
    y[1,t] ~ dpois(x[1,t])
    y[2,t] ~ dpois(x[2,t])
    y[3,t] ~ dpois(x[3,t])
  }
  
  for(t in 1:(N_est-1)){
    # process model
    p[1:3,t] ~ dmnorm(x[1:3,t], SIGMA)
    x[1,t+1] <- max(p[1,t], 0)
    x[2,t+1] <- max(p[2,t], 0)
    x[3,t+1] <- max(p[3,t], 0)
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






