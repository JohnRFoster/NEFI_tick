##   This is the model for Cary Ticks on a Daily scale.    
##   
##   Transition in the questing matrix is estimated with a random intercept.
##
##   Survival in the questing matrix is estimated with the a random intercept, and a 
##   fixed effect on temperature.
##
##   State transitions to capture phenology are lineary with respect to cumulative gdd
##   Cumulative gdd has been centered and scaled by standard deviation
##
##   Using the Green Control site from Cary.         
##                                                                 
##   The state is estimated for every sampling day only, 
##   but demographic params are estimated daily.                       
##
##   Zero Inflated Poisson Data model by life stage. 

#'@param site.run Which site is running? One of "Green Control", "Henry Conrtol", "Tea Control"
#'@param met.variable One of "temp", "precip", "rh"
#'@param n.adapt number of adaption iterations
#'@param n.chains number of chains, 1 for parallel batch jobs
#'@param burnin number of burnin iterations
#'@param thin thinning interval
#'@param n.iter number of iterations post burnin

run_model <- function(site.run,met.variable,n.adapt,n.chains,burnin,thin,n.iter){
  
  start.time <- Sys.time()
  print(start.time)
  
  sites <- c("Green Control","Henry Control","Tea Control")
  
  data <- cary_ticks_met_JAGS(sites)
  
  data <- site_data_met(site = site.run, met.variable = met.variable, data)
  data$R <- diag(1,3,3)
  data$b0 <- as.vector(c(0,0))      ## regression beta means
  data$Vb <- solve(diag(100,2))   ## regression beta precision
  data$gdd <- as.vector(scale(data$gdd))
  data$year.index <- rep(cumsum(table(data$year)), table(data$year))
  
  
  inits <- function(){list(x = data$y,
                           repro.mu = rnorm(1, 1, 0.001),
                           phi.l.mu = rnorm(1, 3.07431, 0.001),
                           phi.n.mu = rnorm(1, 4.15213, 0.001),
                           phi.a.mu = rnorm(1, 5, 0.001),
                           grow.ln.mu = runif(1, -3.07431, 0.001),
                           grow.na.mu = runif(1, -4.15213, 0.001),
                           beta.11 = rnorm(1, 0.11420, 0.001),
                           beta.22 = rnorm(1, -0.05054, 0.001),
                           beta.33 = rnorm(1, 0, 0.001))}
  
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
               "beta.l2n",
               "beta.n2a",
               "beta.a2l",
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
  
  ### Fixed effects priors
  beta.11 ~ dnorm(0.1142,2104.1999)
  beta.22 ~ dnorm(-0.05054,383.5644)
  beta.33 ~ dnorm(0,0.01)
  beta.l2n ~ dmnorm(b0,Vb)
  beta.n2a ~ dmnorm(b0,Vb)
  beta.a2l ~ dmnorm(b0,Vb)
  
  ## prior probablity we observe ticks by life stage
  theta.larvae ~ dunif(0,1)
  theta.nymph ~ dunif(0,1)
  theta.adult ~ dunif(0,1)
  
  for(t in 1:39){
    met[temp.mis[t]] ~ dunif(-30,23)
  }
  
  ### first latent process 
  for(l in 1:3){
    x[l, 1] ~ dpois(2) 
  } # l
  
  ### define parameters
  for(t in 1:N_days){   # loop over every day in time series
  
    ## site specific  daily transtion matrices;
    ## parameters are a linear combination of
    ## the global parameter (intercept) and a 
    ## site random effect
    
    logit(A.day[1,1,t]) <- phi.l.mu + beta.11*met[t]
    logit(A.day[2,1,t]) <- grow.ln.mu
    logit(A.day[2,2,t]) <- phi.n.mu + beta.22*met[t]
    logit(A.day[3,2,t]) <- grow.na.mu
    logit(A.day[3,3,t]) <- phi.a.mu + beta.33*met[t]
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
  
  ### Process Model
  
  for(t in 1:(N_est-1)){
    # molting rates
    logit(l2n[t]) <- beta.l2n[1] + beta.l2n[2]*gdd[t]
    logit(n2a[t]) <- beta.n2a[1] + beta.n2a[2]*gdd[t]
    logit(a2l[t]) <- beta.a2l[1] + beta.a2l[2]*gdd[t]
    
    # number molting
    nymph.molt[t] <- ifelse((gdd[t] >= 500) && (gdd[t] <= 2000), x[1,t] * l2n[t], 0)
    adult.molt[t] <- ifelse((gdd[t] <= 750) || (gdd[t] >= 2500), x[2,t] * n2a[t], 0)
    larva.molt[t] <- ifelse((gdd[t] >= 1500) && (gdd[t] <= 2350), x[3,t] * a2l[t], 0)
    
    # number questing
    Ex[1:3,t] <- TRANS[1:3,1:3,dt.index[t]] %*% x[1:3,t]
    
    # total new ticks
    z[1,t] <- larva.molt[t] + Ex[1,t]
    z[2,t] <- nymph.molt[t] + Ex[2,t]
    z[3,t] <- adult.molt[t] + Ex[3,t]
    
    # process error
    x[1:3,t+1] ~ dmnorm(z[1:3,t], SIGMA)
    
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
  }" 
  
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
  #return < list(out = out, j.model = j.model)
  return(out)
  
  cat("\nJAGS model split and returned\n")
  
  time <- Sys.time() - start.time 
  cat("\nTotal time\n")
  print(time)

}



# 
# print("-------------- DIC ----------------")
# 
# dic.samples(j.model, 10000)




