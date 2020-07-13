
ua_parts <- function(params, ic, type, process.type = "multi"){
  
  Nmc <- nrow(params)
  ua.parts <- list()
  
  if(process.type == "multi"){pattern <- "SIGMA"} # for removing from param matrix
  if(process.type == "independent"){pattern <- "tau"}
    
  # parameter mcmc draws (without process error)
  p.mcmc <- params[,-grep(pattern, colnames(params))]
  
  # parameter names
  params.vec <- colnames(p.mcmc)
  
  for(i in seq_along(params.vec)){
    if("parameter" %in% type){ # with parameter uncertainty
      ua.parts[[i]] <- p.mcmc[,i]
    } else { # without parameter uncertainty
      param.mean <- apply(p.mcmc, 2, mean) # mean parameters
      ua.parts[[i]] <- rep(param.mean[i], Nmc)
    }
  }
  
  # names list parts
  names(ua.parts) <- params.vec
  
  # select appropriate initial conditions
  if("ic" %in% type){
    IC <- ic
  } else {
    IC.mean <- apply(ic, 2, mean)
    IC.names <- names(IC.mean)
    IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
    colnames(IC) <- IC.names
  }
  
  # add IC to ua.parts
  ua.parts$IC <- IC
  
  # process error coloumns 
  if(process.type == "multi"){
    # select appropraite covariance matrix
    if("process" %in% type){ # with process error
      SIGMA <- array(NA, dim = c(3,3,Nmc))
      SIGMA[1,1,] <- params[,"SIGMA[1,1]"]
      SIGMA[1,2,] <- params[,"SIGMA[1,2]"]
      SIGMA[1,3,] <- params[,"SIGMA[1,3]"]
      SIGMA[2,1,] <- params[,"SIGMA[2,1]"]
      SIGMA[2,2,] <- params[,"SIGMA[2,2]"]
      SIGMA[2,3,] <- params[,"SIGMA[2,3]"]
      SIGMA[3,1,] <- params[,"SIGMA[3,1]"]
      SIGMA[3,2,] <- params[,"SIGMA[3,2]"]
      SIGMA[3,3,] <- params[,"SIGMA[3,3]"]
      
      # convert from precision to standard dev
      for(i in 1:Nmc){
        SIGMA[,,i] <- solve(SIGMA[,,i])
      }
    } else { # without process error
      SIGMA <- array(0, dim = c(3,3,Nmc))
    }
    ua.parts$SIGMA <- SIGMA # add sigma to ua.parts
  } 
  
  if(process.type == "independent"){
    if("process" %in% type){ # with process error
      ua.parts$tau.l <- sqrt(1 / params[,"tau_l"])
      ua.parts$tau.n <- sqrt(1 / params[,"tau_n"])
      ua.parts$tau.a <- sqrt(1 / params[,"tau_a"])
    } else {
      ua.parts$tau.l <- rep(0, Nmc)
      ua.parts$tau.n <- rep(0, Nmc)
      ua.parts$tau.a <- rep(0, Nmc)
    }
  }
  
  return(ua.parts)

}
