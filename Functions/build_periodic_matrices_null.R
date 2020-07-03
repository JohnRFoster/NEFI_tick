library(purrr)
library(boot)

build_periodic_matrices_null <- function(ua, data, hb.site = NULL){
  
  Nmc <- length(ua$phi.a.mu)
  
  # pull needed data from list
  N_days <- max(data$N_days)
  gdd <- data$gdd[1:N_days]
  
  if(is.null(hb.site)){
    alpha.a <- rep(0, Nmc)
    alpha.n <- rep(0, Nmc)
    alpha.l <- rep(0, Nmc)
  } else {
    if(hb.site == "Green Control"){s <- 1}
    if(hb.site == "Henry Control"){s <- 2}
    if(hb.site == "Tea Control"){s <- 3}
    if("alpha.a[1]" %in% names(ua)){
      alpha.a <- paste0("alpha.a[", s, "]")
      alpha.a <- pluck(ua, alpha.a)
    } else {
      alpha.a <- rep(0, Nmc)
    }
    if("alpha.n[1]" %in% names(ua)){
      alpha.n <- paste0("alpha.n[", s, "]")
      alpha.n <- pluck(ua, alpha.n)
    } else {
      alpha.n <- rep(0, Nmc)
    }
    if("alpha.l[1]" %in% names(ua)){
      alpha.l <- paste0("alpha.l[", s, "]")
      alpha.l <- pluck(ua, alpha.l)
    } else {
      alpha.l <- rep(0, Nmc)
    }
  }
  
  if("rho.l" %in% names(ua)){
    rho.l <- ua$rho.l
  } else {
    rho.l <- rep(1500, Nmc)
  }
  if("rho.n" %in% names(ua)){
    rho.n <- ua$rho.n
  } else {
    rho.n <- rep(500, Nmc)
  }
  if("rho.a" %in% names(ua)){
    rho.a <- ua$rho.a
  } else {
    rho.a <- rep(2500, Nmc)
  }
    
  # storage
  A <- array(0, dim=c(3,3,N_days,Nmc))
  
  # build matrices
  for(m in seq_along(ua$phi.a.mu)){
    
    phi.11 <- inv.logit(ua$phi.l.mu[m])
    phi.22 <- inv.logit(ua$phi.n.mu[m])
    lambda <- ifelse(gdd >= rho.l[m] & gdd <= 2500, ua$repro.mu[m], 0)
    theta.32 <- ifelse(gdd <= 1000 | gdd >= rho.a[m], inv.logit(ua$grow.na.mu[m]), 0)
    theta.21 <- ifelse(gdd >= rho.n[m] & gdd <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)  
    
    # draw transition matrix with
    
    A[1,1,,m] <- phi.11*(1-theta.21) + alpha.l[m]
    A[2,1,,m] <- phi.11*theta.21
    A[2,2,,m] <- phi.22*(1-theta.32) + alpha.n[m]
    A[3,2,,m] <- phi.22*theta.32
    A[3,3,,m] <- inv.logit(ua$phi.a.mu[m] + alpha.a[m])
    A[1,3,,m] <- lambda
  }
  return(A) 
}