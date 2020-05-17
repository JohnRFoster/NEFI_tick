mouse_forecast <- function(days, params, ic, precip, temp){
  
  library(boot)
  library(mvtnorm)
  
  Nmc <- nrow(params)
  
  # hindcast/forecast settings
  dt <- 1:days

  ind <- nrow(ic)
  
  # storage
  
  # dim = c(latent or observed, ensemble member, individual, days)
  n <-  array(0, dim = c(2, ind, Nmc, days))
  
  lambda.mean <- params[, "lambda.mean"]
  beta.1 <- params[, "beta[1]"]
  beta.2 <- params[, "beta[2]"]
  theta <- params[, "theta"]
  gamma <- params[, grep("gamma[", colnames(params), fixed = TRUE)]
  
  for(m in 1:Nmc){
    # counter
    # if(m%%100==0){cat(round(m/Nmc*100), "% ensembles completed\n")}
    
    lambda <- inv.logit(lambda.mean[m] + beta.1[m]*precip + beta.2[m]*temp)

    phi <- rep(NA, days)
    for(d in 1:days){  ## loop over capture event
      phi[d] <- prod(lambda[1:dt[d]])
    } # d
    
    # initialize condition
    IC <- ic[,m]
    q.IC <- 1 - IC
    n[1,,m,1] <- IC
    
    # for(i in 1:nrow(q)){
    #   q.prod[i] <- prod(q[i, 1:(t-1)])
    # }
    
    for(t in 2:days){
      
      ## State Process
      mu1 <- phi[t-1] * IC + median(gamma[m,]) * q.IC 
      n[1,,m,t] <- rbinom(ind, 1, mu1)  
      
      ## Observation process
      mu2 <- theta[m] * n[1,,m,t]
      n[2,,m,t] <- rbinom(ind, 1, mu2)  
      
    } # t
    
  }
  return(n)
}