library(ecoforecastR)
library(boot)

source("Functions/mouse_data_jags.R")

mouse_one <- function(site, jags.out, type, data){

  Nmc <- nrow(jags.out)
  
  # get data
  days <- data$days
  dt <- data$dt
  ind <- data$ind
  time <- data$time
  precip <- data$precip
  temp <- data$temp
  temp[data$temp.mis] <- mean(temp, na.rm = TRUE)
  
  # draw parameters
  if("parameter" %in% type){
    lambda.mean <- jags.out[, "lambda.mean"]
    beta.1 <- jags.out[, "beta[1]"]
    beta.2 <- jags.out[, "beta[2]"]
    theta <- jags.out[, "theta"]
    gamma <- jags.out[, grep("gamma[", colnames(jags.out), fixed = TRUE)]
  } else {
    lambda.mean <- rep(mean(jags.out[,"lambda.mean"]), Nmc)
    beta.1 <- rep(mean(jags.out[, "beta[1]"]), Nmc)
    beta.2 <- rep(mean(jags.out[, "beta[2]"]), Nmc)
    theta <- rep(mean(jags.out[, "theta"]), Nmc)
    gamma.vec <- colMeans(jags.out[, grep("gamma[", colnames(jags.out), fixed = TRUE)])
    gamma <- rbind(gamma.vec, gamma.vec)
    for(m in 2:Nmc){
      gamma <- rbind(gamma, gamma.vec)
    }
  }
  
  # draw initial conditions
  if("ic" %in% type){
    N.jags <- jags.out[, grep("N[", colnames(jags.out), fixed = TRUE)]  
  } else {
    N.jags.vec <- round(colMeans(jags.out[, grep("N[", colnames(jags.out), fixed = TRUE)]))
    N.jags <- rbind(N.jags.vec, N.jags.vec)
    for(m in 2:Nmc){
      N.jags <- rbind(N.jags, N.jags.vec)
    }
  }
  
  # storage
  y <- p <- array(0, dim = c(ind, time, Nmc))
  q <- matrix(0, ind, time+1)
  latent.state <- state.obs <- matrix(NA, Nmc, time)
  
  for(m in 1:Nmc){
    
    # % counter
    percent <- round(m/Nmc*100)
    if(percent%%20==0){cat(percent, "% ensembles completed\n")}
    
    lambda <- inv.logit(lambda.mean[m] + beta.1[m]*precip[,1] + beta.2[m]*temp[,1])
    
    phi <- rep(NA, time)
    for(d in 1:(time-1)){  ## loop over capture event
      phi[d] <- prod(lambda[dt[d]:(dt[d+1])])
    } # d
    
    for(t in 2:(time-1)){
      
      # initialize
      p[1:N.jags[m,t-1],t-1,m] <- 1
      
      ## State Process
      q[,t-1] <- 1 - p[,t-1,m] 
      q.prod <- rep(NA, nrow(q))
      
      # q.prod <- apply(q, 1, function(x) prod(x[1:(t-1)]))
      
      for(i in 1:nrow(q)){
        q.prod[i] <- prod(q[i, 1:(t-1)])
      }
      
      mu1 <- phi[t-1] * p[,t-1,m] + gamma[m,t] * q.prod
      
      if("process" %in% type){
        # Bernoulli trial
        p[,t,m] <- rbinom(ind, 1, mu1*p[,t-1,m])  
      } else {
        p[,t,m] <- mu1*p[,t-1,m]
      }
      
      ## Observation process
      mu2 <- theta[m] * p[,t,m]
      
      if("observation" %in% type){
        # Bernoulli trial
        y[,t,m] <- rbinom(ind, 1, mu2)  
      } else {
        y[,t,m] <- mu2
      }
    } # t
    
    latent.state[m,] <- apply(p[,,m], 2, sum)
    state.obs[m,] <- apply(y[,,m], 2, sum)
  }
  return(list(latent.state = latent.state,
              latent.state.obs = state.obs))
}





