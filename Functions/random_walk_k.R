random_walk_k <- function(params, ic, N_est, start){
  
  Nmc <- nrow(params)
  
  # covariance matrix
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
  
  pred.seq <- start:(N_est-1)
  pred <- array(NA, dim = c(3, N_est, Nmc))
  for(m in 1:Nmc){
    if(m %% 500 == 0) cat(m, "mcmc samples complete\n")
    for(t in pred.seq){
      
      if(t == start){
        l <- ic[m, paste0("x[1,",start,"]")]
        n <- ic[m, paste0("x[2,",start,"]")]
        a <- ic[m, paste0("x[3,",start,"]")] 
        obs <- c(l, n, a)
        pred[,start,m] <- obs
      } else {
        obs <- pred[,t,m]
      }
      
      # process model
      # p <- rmvnorm(1, obs,solve(SIGMA[,,m])) # convert from precision to covariance
      p <- rmvnp(1, obs, SIGMA[,,m])
      ex.p <- rep(NA, 3)
      ex.p[1] <- max(p[1], 0)
      ex.p[2] <- max(p[2], 0)
      ex.p[3] <- max(p[3], 0)
      
      # observation model
      pred[, t+1, m] <- rpois(3, ex.p)
    }  
  }
  return(pred)
}