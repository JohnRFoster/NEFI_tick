
ua_parts <- function(params, ic){
  ua.parts <- list()
  
  # process error coloumns 
  pattern <- "SIGMA"
  
  # parameter mcmc draws (without process error)
  p.mcmc <- params[,-grep(pattern, colnames(params))]
  
  # mean parameters
  param.mean <- apply(p.mcmc, 2, mean)
  
  # parameter names
  params.vec <- colnames(p.mcmc)
  
  for(i in seq_along(params.vec)){
    if("parameter" %in% type){ # with parameter uncertainty
      ua.parts[[i]] <- p.mcmc[,i]
    } else { # without parameter uncertainty
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
   
  # add sigma to ua.parts
  ua.parts$SIGMA <- SIGMA
  
  return(ua.parts)

}
