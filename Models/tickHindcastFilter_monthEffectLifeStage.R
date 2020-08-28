

# for updating priors, need to extend for other known.vals
update_data <- function(params.mat, month.seq, known.vals = NULL){
  
  # moment match beta distribution
  beta_moment_match <- function(param){
    ex <- mean(param)
    v <- var(param)
    alpha <- (((ex^2) * (1-ex)) / v) - ex
    beta <- (((ex*(1-ex)) / v) - 1) * (1-ex)
    moments <- list(alpha = alpha, beta = beta)
    return(moments)
  }
  
  if(is.null(known.vals)){
    data <- list()
  } else {
    data <- known.vals  
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
  
  data$month.mu <- matrix(0, 3, 12) 
  data$month.prec <- matrix(0.001, 3, 12)
  data$n.months <- 1:12
  for(month in month.seq){
    for(s in 1:3){
      month.index <- paste0("alpha.month[", s, ",", month, "]")
      data$month.mu[s, month] <- mean(params.mat[,month.index])   
      data$month.prec[s, month] <- 1 / var(params.mat[,month.index])  
    }
  }
  
  if(!("SIGMA" %in% names(known.vals))){
    wish_df <- function(Om, X, i, j, col) {
      (Om[i, j]^2 + Om[i, i] * Om[j, j]) / var(X[, col])
    }
    
    # process error draws, precision
    OMEGA <- params.mat[,grep("SIGMA", colnames(params.mat))]
    
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
  } 
  return(data)
}

#'@param site.run Which site is running? One of "Green Control", "Henry Conrtol", "Tea Control"
#'@param n.adapt number of adaption iterations
#'@param n.chains number of chains, 1 for parallel batch jobs
#'@param burnin number of burnin iterations
#'@param thin thinning interval
#'@param n.iter number of iterations post burnin

run_jagsFilter <- function(data, n.adapt, n.chains,
                           known.vals = NULL) {
  
  
  inits <- function() {
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
  
  model = " model {

  ### global priors
  phi.l.mu ~ dnorm(larva.mu, larva.prec)               # larvae survival
  phi.n.mu ~ dnorm(nymph.mu, nymph.prec)               # nymph survival
  phi.a.mu ~ dnorm(adult.mu, adult.prec)               # adult survival
  grow.ln.mu ~ dnorm(l2n.mu, l2n.prec)             # larvae -> nymph transition
  grow.na.mu ~ dnorm(n2a.mu, n2a.prec)             # nymph -> adult transition
  repro.mu ~ dnorm(a2l.mu, a2l.prec) T(0,)              # adult -> larvae transition (reproduction)
  rho.l ~ dnorm(rho.l.mu, rho.l.prec) T(0,)
  rho.n ~ dnorm(rho.n.mu, rho.n.prec) T(0,)
  rho.a ~ dnorm(rho.a.mu, rho.a.prec) T(0,)
  
  ### precision priors
  SIGMA ~ dwish(R, k)         # mvn [3 x 3] site process
  
  ### random month prior
  for(s in 1:3){
    for(month in n.months){
      alpha.month[s,month] ~ dnorm(month.mu[s, month], month.prec[s, month])
    }
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
  for(t in 1:(n.days)){
  
    # expected number questing
    Ex[1:3,t] <- TRANS[1:3,1:3,t] %*% x[1:3,t] + alpha.month[1:3, month.index[t]]
  
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
  
  
  load.module("glm")
  
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


