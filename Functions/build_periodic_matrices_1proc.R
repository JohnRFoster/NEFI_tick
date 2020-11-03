build_periodic_matrices_1proc <- function(ua, data, mice.site = NA){
  
  # pull needed data from list
  N_days <- data$N_days
  gdd <- data$gdd[1:N_days]
  
  Nmc <- length(ua$phi.a.mu)
  
  # storage
  A <- array(0, dim=c(3,3,N_days,Nmc))
  
  # get mice data
  if(!is.na(mice.site)){
    source("/projectnb/dietzelab/fosterj/NEFI_tick/Functions/mice_estimated_jags.R")
    cat("Getting mice from", mice.site, "\n")
    mice <- mice_estimated_jags(mice.site)
    mice.mean <- mice$mice.mean
    mice.prec <- mice$mice.prec
  }
  
  # build matrices
  for(m in seq_along(ua$phi.a.mu)){
    
    if("beta.l" %in% names(ua)){
      met <- data$met[1:N_days]
      met[data$met.mis] <- mean(met, na.rm = TRUE)
      phi.11 <- boot::inv.logit(ua$phi.l.mu[m] + ua$beta.l[m]*met)
    } else {
      phi.11 <- boot::inv.logit(ua$phi.l.mu[m])
    }
    
    if("beta.n" %in% names(ua)){
      met <- data$met[1:N_days]
      met[data$met.mis] <- mean(met, na.rm = TRUE)
      phi.22 <- boot::inv.logit(ua$phi.n.mu[m] + ua$beta.n[m]*met)
    } else {
      phi.22 <- boot::inv.logit(ua$phi.n.mu[m])
    }
    
    if("beta.a" %in% names(ua)){
      met <- data$met[1:N_days]
      met[data$met.mis] <- mean(met, na.rm = TRUE)
      phi.33 <- boot::inv.logit(ua$phi.a.mu[m] + ua$beta.a[m]*met)
    } else {
      phi.33 <- boot::inv.logit(ua$phi.a.mu[m])
    }
    
    if(is.na(mice.site)){
      l2n <- boot::inv.logit(ua$grow.ln.mu[m])
      n2a <- boot::inv.logit(ua$grow.na.mu[m])
    } else {  
      mice.real <- rep(NA, length(mice.mean))
      for(mm in seq_along(mice.mean)){
        mice.real[mm] <- rnorm(1, mice.mean[mm], LaplacesDemon::prec2sd(mice.prec[mm]))
      }
      if("beta.l2n" %in% names(ua)){
        l2n <- boot::inv.logit(ua$grow.ln.mu[m] + ua$beta.l2n[m]*mice.real)
      } else {
        l2n <- boot::inv.logit(ua$grow.ln.mu[m])
      }
      if("beta.n2a" %in% names(ua)){
        n2a <- boot::inv.logit(ua$grow.na.mu[m] + ua$beta.n2a[m]*mice.real)
      } else {
        n2a <- boot::inv.logit(ua$grow.na.mu[m])
      }
    }
      
    lambda <- ifelse(gdd >= 1500 & gdd <= 2500, ua$repro.mu[m], 0)
    theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, n2a, 0)
    theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, l2n, 0)  
    
    # draw transition matrix
    A[1,1,,m] <- phi.11*(1-theta.21)
    A[2,1,,m] <- phi.11*theta.21
    A[2,2,,m] <- phi.22*(1-theta.32)
    A[3,2,,m] <- phi.22*theta.32
    A[3,3,,m] <- phi.33
    A[1,3,,m] <- lambda
  }
  return(A)
}
