##   This is the one step ahead prediction for temp driving survival model  
##
##   Uncertainty partitioned for each parameter as well as initial condition and process     
##                                                                 
##   The state is estimated for every sampling day only, 
##   but demographic params are estimated daily.                       
##
##   Zero Inflated Poisson Data model by life stage. 

#' @param type type of uncertainty to simulate. Can be any of "deterministic","ic","parameter","process","larvae survival","nymph survival","adult survival","larvae-to-nymph","nymph-to-adult","reproduction"
#' @param site site, one of "Green Control","Henry Control","Tea Control"
#' @param params matrix of JAGS output of parameters
#' @param ic matrix of JAGS output of states
#' @param data model data fed to JAGS
#' @param Nmc Number of mc samples
#' @param draw Vector of random draws (row numbers) or sampling
#' @export    

predict_state_one_temp_matrix <- function(type, site, params, ic, data, Nmc, draw){
  
  # call data from model run for storage dimensions and indexing
  if(site == "Green Control"){
    s <- 1
  } else if(site == "Henry Control"){
    s <- 2
  } else {
    s <- 3
  }
  N_est <- data$N_est[s]
  N_days <- data$N_days[s]
  dt.index <- data$dt.index[s,]
  df <- data$df[s,]
  met <- data$met[1:N_days,,s]
  temp.mis <- data$temp.mis[,s]
  rh.mis <- data$rh.mis[,s]
  
  temp <- met[,1]
  mean.temp <- mean(temp, na.rm = TRUE)
  for(i in 1:length(temp.mis)){
    temp[temp.mis[i]] <- mean.temp
  }
  met.x <- temp

  # storage
  pred <- array(dim = c(3,N_est,Nmc))
  A <- array(0, dim = c(3,3,N_days))
  
  # mean parameters
  param.mean <- apply(params, 2, mean)
  
  IC.mean <- apply(ic, 2, mean)
  IC.names <- names(IC.mean)
  IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
  colnames(IC) <- IC.names
  
  alpha.11 <- params[,paste("alpha.11[", s, "]", sep = "")]
  alpha.13 <- params[,paste("alpha.13[", s, "]", sep = "")]
  alpha.22 <- params[,paste("alpha.22[", s, "]", sep = "")]
  alpha.33 <- params[,paste("alpha.33[", s, "]", sep = "")]
  alpha.32 <- params[,paste("alpha.32[", s, "]", sep = "")]
  
  # select appropriate matrix parameters
  if("larvae survival" %in% type){
    phi.l.mu <- params[draw,"phi.l.mu"]
    beta.11 <- params[draw, "beta.11"]
    alpha.11 <- alpha.11[draw]
  } else {
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    beta.11 <- rep(param.mean["beta.11"], Nmc)
    alpha.11 <- rep(mean(alpha.11), Nmc)
  }
  
  if("nymph survival" %in% type){
    phi.n.mu <- params[draw,"phi.n.mu"]
    beta.22 <- params[draw, "beta.22"]
    alpha.22 <- alpha.22[draw]
  } else {
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    beta.22 <- rep(param.mean["beta.22"], Nmc)
    alpha.22 <- rep(mean(alpha.22), Nmc)
  }
  
  if("adult survival" %in% type){
    phi.a.mu <- params[draw,"phi.a.mu"]
    beta.33 <- params[draw, "beta.33"]  
    alpha.33 <- alpha.33[draw]
  } else {
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    beta.33 <- rep(param.mean["beta.33"], Nmc)
    alpha.33 <- rep(mean(alpha.33), Nmc)
  }
  
  if("larva-to-nymph" %in% type){
    grow.ln.mu <- params[draw,"grow.ln.mu"]
  } else {
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
  }
  
  if("nymph-to-adult" %in% type){
    grow.na.mu <- params[draw,"grow.na.mu"]
    alpha.32 <- alpha.32[draw]
  } else {
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    alpha.32 <- rep(mean(alpha.32), Nmc)
  }
    
  if("reproduction" %in% type){
    repro.mu <- params[draw,"repro.mu"]
    alpha.13 <- alpha.13[draw]
  } else {
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
    alpha.13 <- rep(mean(alpha.13), Nmc)
  }
  
  SIGMA <- array(0, dim = c(3,3,Nmc))

  # run mcmc sampling
  for(m in 1:Nmc){
    
    # draw transition matrix
    A[1,1,] <- inv.logit(phi.l.mu[m] + beta.11[m]*met.x + alpha.11[m])
    A[1,3,] <- exp(repro.mu[m] + alpha.13[m])
    A[2,1,] <- inv.logit(grow.ln.mu[m])
    A[2,2,] <- inv.logit(phi.n.mu[m] + beta.22[m]*met.x + alpha.22[m])
    A[3,2,] <- inv.logit(grow.na.mu[m] + alpha.32[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + beta.33[m]*met.x + alpha.33[m])
    
    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){                 # loop over the number of sampling days - 1
      if(t == 1){
        TRANS <- A[,,1] %*% A[,,2]
        for(d in 3:df[1]){
          TRANS <- TRANS %*% A[,,d]
        }
      } else {
        for(d in 1:df[t]){
          if(d == 1){
            TRANS <- A[,,dt.index[t-1]] %*% A[,,dt.index[t-1]+1]
          } else {
            TRANS <- TRANS %*% A[,,dt.index[t-1]+d]
          }
        }
      }
      
      ## initial condition
      l <- IC[m, paste("x[1,",t, ",", s,"]",sep="")]
      n <- IC[m, paste("x[2,",t, ",", s,"]",sep="")]
      a <- IC[m, paste("x[3,",t, ",", s,"]",sep="")]
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
      pred[1,t,m] <- max(est.mvnorm[1], 0)
      pred[2,t,m] <- max(est.mvnorm[2], 0)
      pred[3,t,m] <- max(est.mvnorm[3], 0)
    }
  }
  return(pred)
} # close function
