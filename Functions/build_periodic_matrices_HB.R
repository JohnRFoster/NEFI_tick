library(boot)

build_periodic_matrices_HB <- function(ua, data){
  
  Nmc <- length(ua$repro.mu)
  A.list <- list() # storage
  
  for(s in 1:3){
    
    met <- data$met
    met[data$met.mis] <- mean(met, na.rm = TRUE)
    
    day.seq <- data$dt.index[s,1]:data$dt.index[s,data$N_est[s]]
    gdd <- data$gdd[day.seq]
    met <- met[day.seq]
    
    # adult survival random effect, ran model two ways so need grab from ua
    # depending on how the bayes model was coded
    if("alpha.a[1]" %in% names(ua)){
      alpha <- ua[[paste0("alpha.a[", s, "]")]]
    } else {
      alpha <- rep(0, Nmc)
    }
    
    if("phi.a.mu[1]" %in% names(ua)){
      phi.a.mu <- ua[[paste0("phi.a.mu[", s, "]")]]
    } else {
      phi.a.mu <- ua$phi.a.mu
    }
    
    # build matrices
    A <- array(0, dim=c(3,3,length(day.seq),Nmc)) # storage
    for(m in seq_along(ua$repro.mu)){
      
      if("beta.l" %in% names(ua)){
        phi.11 <- inv.logit(ua$phi.l.mu[m] + ua$beta.l[m]*met)
      } else {
        phi.11 <- inv.logit(ua$phi.l.mu[m])
      }
      
      if("beta.n" %in% names(ua)){
        phi.22 <- inv.logit(ua$phi.n.mu[m] + ua$beta.n[m]*met)
      } else {
        phi.22 <- inv.logit(ua$phi.n.mu[m])
      }
      
      if("beta.a" %in% names(ua)){
        phi.33 <- inv.logit(phi.a.mu[m] + ua$beta.a[m]*met + alpha[m])
      } else {
        phi.33 <- inv.logit(phi.a.mu[m] + alpha[m])
      }
      
      lambda <- ifelse(gdd >= 1500 & gdd <= 2500, ua$repro.mu[m], 0)
      theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, inv.logit(ua$grow.na.mu[m]), 0)
      theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)  
      
      # draw transition matrix
      A[1,1,,m] <- phi.11*(1-theta.21)
      A[2,1,,m] <- phi.11*theta.21
      A[2,2,,m] <- phi.22*(1-theta.32)
      A[3,2,,m] <- phi.22*theta.32
      A[3,3,,m] <- phi.33 
      A[1,3,,m] <- lambda
    }  
    A.list[[s]] <- A
  }
  return(A.list)
}
