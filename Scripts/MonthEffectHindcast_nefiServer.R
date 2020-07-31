# =================================================== #
#   This script runs a hindcast for a single site at  #
#   the Cary Institute of Ecosystem Studies           #
#                                                     #
#   The method for data assimilation is a jags filter #
#                                                     #
#   This model has a fixed month effect that is the   #
#   same for all life stages, no met on survival      #
# =================================================== #

library(rjags)
library(tidyverse)
library(lubridate)

# =================================================== #
#               Data and JAGS inputs                  #
# =================================================== #

n.adapt <- 2000       # adaptation iterations
n.iter <- 10000       # sampling iterations
n.chains <- 3         # number of chains
iter2save <- 10000    # number of samples to save (thinning)

# grid <- "Henry Control" # site name
# met.driver <- NULL      # this model does not use met to drive survival

# load and parse data 
data.orig <- read.csv("DataHindcast_HenryControl.csv")

# met for everyday during the hindcast period
met <- data.orig[complete.cases(data.orig[,2:4]), 2:4]
met$date <- ymd(met$date)

# tick observations and dates of observation
ticks.observed <- t(data.orig[complete.cases(data.orig[,5:7]), 5:7])
days <- data.orig$days.csv[complete.cases(data.orig$days.csv)]
days <- ymd(days)

# mcmc samples from training fit
params.mcmc <- data.orig[,9:ncol(data.orig)]

colnames(params.mcmc) <- 
c("SIGMA[1,1]", # hardcoding these because write.csv doesn't like brackets
  "SIGMA[2,1]",
  "SIGMA[3,1]",
  "SIGMA[1,2]",
  "SIGMA[2,2]",
  "SIGMA[3,2]",
  "SIGMA[1,3]",
  "SIGMA[2,3]",
  "SIGMA[3,3]",
  "alpha.month[4]",
  "alpha.month[5]",
  "alpha.month[6]",
  "alpha.month[7]",
  "alpha.month[8]",
  "alpha.month[9]",
  "alpha.month[10]",
  "alpha.month[11]",
  "alpha.month[12]",
  "beta.l.obs",
  "deviance",
  "grow.ln.mu",
  "grow.na.mu",
  "phi.a.mu",
  "phi.l.mu",
  "phi.n.mu",
  "repro.mu",
  "rho.a",
  "rho.l",
  "rho.n",
  "theta.adult",
  "theta.nymph"
)

# =================================================== #
#               Function and constants                #
# =================================================== #

update_data <- function(params.mat, month.seq){
  # for updating parameter priors after each hindcast
  # params.mat = mcmc samples, matrix of parameters [samples, parameter]
  # month.seq = vector of months - integer
  
  data <- list() # data list for jags
  
  # moment match beta distribution
  beta_moment_match <- function(param){
    ex <- mean(param)
    v <- var(param)
    alpha <- (((ex^2) * (1-ex)) / v) - ex
    beta <- (((ex*(1-ex)) / v) - 1) * (1-ex)
    moments <- list(alpha = alpha, beta = beta)
    return(moments)
  }

  prior.theta.adult <- beta_moment_match(params.mat[,"theta.adult"])
  prior.theta.nymph <- beta_moment_match(params.mat[,"theta.nymph"])
  data$beta.adult.alpha <- prior.theta.adult$alpha
  data$beta.adult.beta <- prior.theta.adult$beta
  data$beta.nymph.alpha <- prior.theta.nymph$alpha
  data$beta.nymph.beta <- prior.theta.nymph$beta
  data$beta.l.obs.mu <- mean(params.mat[,"beta.l.obs"])
  data$beta.l.obs.prec <- 1 / var(params.mat[,"beta.l.obs"])
  data$larva.mu <- mean(params.mat[,"phi.l.mu"])
  data$larva.prec <- 1 / var(params.mat[,"phi.l.mu"])
  data$nymph.mu <- mean(params.mat[,"phi.n.mu"])
  data$nymph.prec <- 1 / var(params.mat[,"phi.n.mu"])
  data$adult.mu <- mean(params.mat[,"phi.a.mu"])
  data$adult.prec <- 1 / var(params.mat[,"phi.a.mu"])
  data$l2n.mu <- mean(params.mat[,"grow.ln.mu"])
  data$l2n.prec <- 1 / var(params.mat[,"grow.ln.mu"])
  data$n2a.mu <- mean(params.mat[,"grow.na.mu"])
  data$n2a.prec <- 1 / var(params.mat[,"grow.na.mu"])
  data$a2l.mu <- mean(params.mat[,"repro.mu"])
  data$a2l.prec <- 1 / var(params.mat[,"repro.mu"])
  data$rho.l.mu <- mean(params.mat[,"rho.l"])
  data$rho.l.prec <- 1 / var(params.mat[,"rho.l"])
  data$rho.n.mu <- mean(params.mat[,"rho.n"])
  data$rho.n.prec <- 1 / var(params.mat[,"rho.n"])
  data$rho.a.mu <- mean(params.mat[,"rho.a"])
  data$rho.a.prec <- 1 / var(params.mat[,"rho.a"])  
  
  # start with flat priors for month effects, then update
  data$month.mu <- rep(0, 12) 
  data$month.prec <- rep(1/10^2, 12)
  data$n.months <- 1:12
  for(month in month.seq){
    month.index <- paste0("alpha.month[", month, "]")
    data$month.mu[month] <- mean(params.mat[,month.index])   
    data$month.prec[month] <- 1 / var(params.mat[,month.index])
  }
  
  # update wishart
  wish_df <- function(Om, X, i, j, col) {
    (Om[i, j]^2 + Om[i, i] * Om[j, j]) / var(X[, col])
  }
  
  OMEGA <- params.mat[,grep("SIGMA", colnames(params.mat))] # process error draws, precision
  q.bar <- matrix(apply(OMEGA, 2, mean), 3, 3)  # Mean Omega, Precision
  
  col <- matrix(1:3 ^ 2, 3, 3)
  WV <- matrix(NA, 3, 3)
  for(r in 1:3){
    for(c in 1:3){
      WV[r, c] <- wish_df(q.bar, OMEGA, r, c, col[r, c])
    }
  }
  
  data$R <- solve(q.bar) * mean(WV)
  data$k <- mean(WV)  
  
  return(data)
}

convergence_check <- function(jags.out, model, monitor,
                              n.iter = 10000, min.eff.size = 5000, GBR.thresh = 1.02){
  # checks converegence after coda.samples, will resample if not converged or
  # sample sizes not big enough
  # jags.out = samples from coda.samples
  # model = compiled model object from jags.model
  # monitor = variables to monitor
  
  enough.samples <- converge <- FALSE
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  
  # convergence check on parameters
  cat("Calculating PSRF\n")
  GBR.vals <- gelman.diag(out$params, multivariate = FALSE)
  GBR.vals
  converge <- max(GBR.vals$psrf) < GBR.thresh
  cat("Convergence:", converge, "\n")
  
  if(converge){
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin)){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      enough.samples <- min(effectiveSize(out$params)) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
    }
  }
  
  counter <- 1
  while(!converge & !enough.samples){
    counter <- counter + 1
    cat("coda.samples call number:", counter, "\n")
    jags.out <- coda.samples(model = model$j.model,
                             variable.names = monitor,
                             n.iter = 10000)
    
    cat("coda samples done, checking mcmc \n")
    
    ## split output
    out <- list(params = NULL, predict = NULL)
    mfit <- as.matrix(jags.out, chains = TRUE)
    pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
    chain.col <- which(colnames(mfit) == "CHAIN")
    out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
    out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
    
    # convergence check on parameters
    cat("Calculating PSRF\n")
    GBR.vals <- gelman.diag(out$params, multivariate = FALSE)
    GBR.vals
    converge <- max(GBR.vals$psrf) < GBR.thresh
    cat("Convergence:", converge, "\n")
    if(!converge) next
    
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin)){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      enough.samples <- min(effectiveSize(out$params)) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
    }
  }
  return(out)
}


# variables to monitor
monitor <- c(
  "x",
  "phi.l.mu",
  "phi.n.mu",
  "phi.a.mu",
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

# define jags model
model.filter <- " model {

  ### global priors
  phi.l.mu ~ dnorm(larva.mu, larva.prec)       # larvae survival
  phi.n.mu ~ dnorm(nymph.mu, nymph.prec)       # nymph survival
  phi.a.mu ~ dnorm(adult.mu, adult.prec)       # adult survival
  grow.ln.mu ~ dnorm(l2n.mu, l2n.prec)         # larvae -> nymph transition
  grow.na.mu ~ dnorm(n2a.mu, n2a.prec)         # nymph -> adult transition
  repro.mu ~ dnorm(a2l.mu, a2l.prec) T(0,)     # adult -> larvae transition (reproduction)
  rho.l ~ dnorm(rho.l.mu, rho.l.prec) T(0,)    # larva phenology start
  rho.n ~ dnorm(rho.n.mu, rho.n.prec) T(0,)    # nymph phenology start
  rho.a ~ dnorm(rho.a.mu, rho.a.prec) T(0,)    # adult phenology start
  
  ### precision priors
  SIGMA ~ dwish(R, k)         # mvn [3 x 3] site process
  
  ### random month prior
  for(month in n.months){
    alpha.month[month] ~ dnorm(month.mu[month], month.prec[month])
  }
  
  ## observation regression priors
  beta.l.obs ~ dnorm(beta.l.obs.mu, beta.l.obs.prec) T(1E-10,)
  theta.nymph ~ dbeta(beta.adult.alpha, beta.adult.beta)
  theta.adult ~ dbeta(beta.nymph.alpha, beta.nymph.beta)
  
  ### first latent process
  x[1, 1] ~ dpois(l.ic)
  x[2, 1] ~ dpois(n.ic)
  x[3, 1] ~ dpois(a.ic)
  
  logit(phi.11) <- phi.l.mu
  logit(phi.22) <- phi.n.mu
  logit(l2n) <- grow.ln.mu
  logit(n2a) <- grow.na.mu
  
  ### define parameters
  for(t in 1:n.days){   # loop over every day in time series
  
    theta.21[t] <- ifelse((gdd[t,1] >= rho.n) && (gdd[t,1] <= 2500),l2n,0)
    theta.32[t] <- ifelse((gdd[t,1] <= 1000) || (gdd[t,1] >= rho.a),n2a,0)
    lambda[t] <- ifelse((gdd[t,1] >= rho.l) && (gdd[t,1] <= 2500),repro.mu,0)
  
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
  TRANS[1:3,1:3,n.days] <- A.day[,,n.days] %*% A.day[,,n.days-1]
  
  for(day in seq.days){
    TRANS[1:3,1:3,day] <- TRANS[1:3,1:3,day+1] %*% A.day[,,day]
  }
  
  ### Process Model
  for(t in 1:(n.days-1)){
  
    # expected number questing
    Ex[1:3,t] <- (TRANS[1:3,1:3,seq.days[t]] %*% x[1:3,t]) + alpha.month[month.index[t]]
  
    # process error
    p[1:3,t] ~ dmnorm(Ex[1:3,t], SIGMA)
    x[1,t+1] <- max(p[1,t], 0)
    x[2,t+1] <- max(p[2,t], 0)
    x[3,t+1] <- max(p[3,t], 0)
  }
  
  ### Data Model ###
  for(t in 1:n.days){
  
    ## fit the blended model to observed data
    y[1,t] ~ dpois(m[1,t])
    y[2,t] ~ dpois(m[2,t])
    y[3,t] ~ dpois(m[3,t])
    
    ## blend the poisson and zero inflation models
    m[1,t] <- x[1,t]*b.larva[t] + 1E-10
    m[2,t] <- x[2,t]*b.nymph[t] + 1E-10
    m[3,t] <- x[3,t]*b.adult[t] + 1E-10
    
    ## observation probability based on temperature
    theta.larva[t] <- 1 / (1 + beta.l.obs*(met.obs[t,1])^2)
    
    ## binary outcome of observation by life stage
    b.larva[t] ~ dbern(theta.larva[t])
    b.nymph[t] ~ dbern(theta.nymph)
    b.adult[t] ~ dbern(theta.adult)
  } # t
}" 

inits <- function() {
  # inits for jags.model
  list(
    repro.mu = max(0, rnorm(1, data$a2l.mu, 1)),
    phi.l.mu = max(0, rnorm(1, data$larva.mu, 1)),
    phi.n.mu = max(0, rnorm(1, data$nymph.mu, 1)),
    phi.a.mu = max(0, rnorm(1, data$adult.mu, 1)),
    rho.a = max(0, rnorm(1, data$rho.a.mu, 100)),
    rho.l = max(0, rnorm(1, data$rho.l.mu, 100)),
    rho.n = max(0, rnorm(1, data$rho.n.mu, 100)),
    grow.ln.mu = min(0, rnorm(1, data$l2n.mu, 1)),
    grow.na.mu = min(0, rnorm(1, data$n2a.mu, 1))
  )
}

# =================================================== #
#                    Run hindcast                     #
# =================================================== #

# data and initial conditions for first hindcast
data <- update_data(params.mcmc, 4:12) # only update month effect params from training fit
data$l.ic <- round(rpois(1, ticks.observed[1,1])) + 1
data$n.ic <- round(rpois(1, ticks.observed[2,1])) + 1
data$a.ic <- round(rpois(1, ticks.observed[3,1])) + 1

hindcast.seq <- 1:(length(days)-1) # for entire time series

# hindcast.seq <- 1:5 # for testing
# t <- 1  # for testing

for (t in seq_along(hindcast.seq)) {
  
  cat("==================================\n")
  cat(round(t/length(hindcast.seq)*100, 2), "% done\n")
  
  ### pull out days ###
  forecast.start.day <- as.character(days[t])  # date forecast issued
  forecast.end.day <- as.character(days[t+1])  # next observation date
  check.day <- as.integer(days[t+1] - days[t]) # day we evaluate forecast and run DA
  
  all.days <- seq.Date(ymd(forecast.start.day)+1, ymd(forecast.end.day), 1)
  data$month.index <- month(all.days) # months in hindcast period
  
  # grab met
  met.subset <- met %>% 
    filter(date > forecast.start.day) %>% 
    filter(date <= forecast.end.day) 
  obs.temp <- met.subset %>% 
    select(min.temp.scale)
  gdd <- met.subset %>% 
    select(cum.gdd)
  
  data$n.days <- nrow(gdd) # number of days in forecast
  data$gdd <- gdd
  data$met.obs <- obs.temp
  
  y <- matrix(NA, 3, check.day)
  y[,1] <- ticks.observed[,t] # only known data is first day of hindcast
  data$y <- y
  data$seq.days <- (data$n.days-1):1
  
  # run filter
  load.module("glm")
  jags.filter <- jags.model(
    file = textConnection(model.filter),
    data = data,
    inits = inits,
    n.adapt = n.adapt,
    n.chains = n.chains
  )
  
  jags.out <- coda.samples(
    model = jags.filter,
    variable.names = monitor,
    n.iter = n.iter
  )
  
  cat("coda samples done, checking mcmc \n")
  
  ## convergence check
  out <- convergence_check(jags.out, jags.filter, monitor)
  
  # =====  this section is just for printing hindcast predictions ===== #
  preds <- as.matrix(out$predict)
  larva <- preds[,seq(1, by = 3, length.out = ncol(preds)/3)]
  nymph <- preds[,seq(2, by = 3, length.out = ncol(preds)/3)]
  adult <- preds[,seq(3, by = 3, length.out = ncol(preds)/3)]
  
  # quantiles for cat statements
  l.q <- round(quantile(larva[,ncol(larva)], c(0.025, 0.5, 0.975)))
  n.q <- round(quantile(nymph[,ncol(nymph)], c(0.025, 0.5, 0.975)))
  a.q <- round(quantile(adult[,ncol(adult)], c(0.025, 0.5, 0.975)))
  
  # cat statements for batch jobs
  cat("Larva oberved:", ticks.observed[1,t+1], "\n")
  cat("Larva forecast CI:", l.q, "\n\n")
  cat("Nymph oberved:", ticks.observed[2,t+1], "\n")
  cat("Nymph forecast CI:", n.q, "\n\n")
  cat("Adult oberved:", ticks.observed[3,t+1], "\n")
  cat("Adult forecast CI:", a.q, "\n\n")
  # =================================================================== #
  
  parameters <- as.matrix(out$params)
  data <- update_data(parameters, 1:12) # update with full posterior, need all months 
  data$l.ic <- round(mean(larva[,ncol(larva)])) # new initial conditions
  data$n.ic <- round(mean(nymph[,ncol(nymph)]))
  data$a.ic <- round(mean(adult[,ncol(adult)]))
  
  # thin for saving
  thin <- seq(1, nrow(parameters), length.out = iter2save)
  preds <- preds[thin,]            # posterior state samples
  parameters <- parameters[thin,]  # posterior param samples
  
  # save as .RData
  # outname <- paste(t, "Tick_hindcast_jagsFilter.RData", sep = "_") # name file
  # save(preds, parameters, data,
  #      file = file.path(out.dir, outname))
  
}

  
cat("=== END ===\n")









