sensitivity_elasticity <- function(A, data, Nmc){
  
  N_est <- data$N_est
  N_days <- data$N_day
  dt.index <- data$dt.index
  df <- data$df
  gdd <- data$gdd[1:N_days]
  
  seq.days <- matrix(NA, N_est-1, max(df, na.rm = TRUE))
  for(i in 1:(N_est-1)){
    xx <- (dt.index[i+1]-1):dt.index[i]
    seq.days[i,1:length(xx)] <- xx
  }
  
  # Product of phase matrices A excluding matrix at time 0
  D <- array(0, dim = c(3,3,N_est-1,Nmc))
  A.agg <- array(0, dim = c(3,3,N_days,Nmc))
  
  for(m in 1:Nmc){
    for(t in 1:(N_est-1)){   # loop over the number of sampling days - 1
      for(day in seq.days[t,1]){
        A.agg[1:3,1:3,day,m] <- A[,,day,m] %*% A[,,day-1,m]
      } 
      for(day in seq.days[t,2:(df[t]-1)]){
        A.agg[1:3,1:3,day,m] <- A.agg[1:3,1:3,day+1,m] %*% A[,,day,m]
      }
      D[1:3,1:3,t,m] <- A.agg[1:3,1:3,day,m]
    }
  }
  
  # perdioc matrix at sampling events; time 0
  Bh <- A[,,dt.index[1:(N_est-1)],]
  
  # storage for sensitivity and elasticity matrices
  S.Bh <- E.Bh <- array(0, dim = dim(D))
  for(m in 1:Nmc){
    for(t in 1:72){
      
      # A.h is the product of all phase matrices within time frame
      A.h <- D[,,t,m] %*% Bh[,,t,m]
      
      # sensitivity of A.h
      S.Ah <- sensitivity(A.h)
      
      # sensitivity of Bh
      S.Bh[,,t,m] <- t(D[,,t,m]) %*% S.Ah
      
      # elasticity of Bh
      E.Bh[,,t,m] <- (1/lambda(Bh[,,t,m]))*Bh[,,t,m]*S.Bh[,,t,m]
    }  
  }
  
  return(list(sensitivity = S.Bh,
              elasticity = E.Bh))
  
}
