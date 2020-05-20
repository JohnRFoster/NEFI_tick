mouse_forecast_parallel <- function(Nmc.seq){
  
  library(boot)
  library(mvtnorm)
  
  Nmc <- length(Nmc.seq)
  # Nmc.per.node <- Nmc / n.slots
  # Nmc.seq <- ((Nmc.per.node * x) - Nmc.per.node + 1):(Nmc.per.node * x)
  
  # params <- params[Nmc.seq,]
  
  # hindcast/forecast settings
  ind <- nrow(ic)
  dt <- 1:n.days
  
  # storage
  
  # dim = c(latent or observed, ensemble member, individual, days)
  n <-  array(0, dim = c(2, ind, Nmc, n.days))
  
  # reduce to Nmc samples to run
  params <- params[Nmc.seq,]
  ic <- ic[,Nmc.seq]
  
  lambda.mean <- params[, "lambda.mean"]
  beta.1 <- params[, "beta[1]"]
  beta.2 <- params[, "beta[2]"]
  theta <- params[, "theta"]
  gamma <- params[, grep("gamma[", colnames(params), fixed = TRUE)]
  
  for(m in 1:Nmc){
    # counter
    # if(m%%100==0){cat(round(m/Nmc*100), "% ensembles completed\n")}
    
    lambda <- inv.logit(lambda.mean[m] + beta.1[m]*precip[,1] + beta.2[m]*temp[,1])

    phi <- rep(NA, n.days)
    for(d in 1:n.days){  ## loop over capture event
      phi[d] <- prod(lambda[1:dt[d]])
    } # d
    
    # initialize condition
    IC <- ic[,m]
    q.IC <- 1 - IC
    n[1,,m,1] <- IC
    n[2,,m,1] <- IC
    
    # for(i in 1:nrow(q)){
    #   q.prod[i] <- prod(q[i, 1:(t-1)])
    # }
    
    for(t in 2:n.days){
      
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