build_periodic_matrices_1proc <- function(ua, data){
  
  # pull needed data from list
  N_days <- data$N_days
  gdd <- data$gdd[1:N_days]
  met <- data$met[1:N_days]
  met[data$met.mis] <- mean(met, na.rm = TRUE)

  Nmc <- length(ua$phi.a.mu)
  
  # storage
  A <- array(0, dim=c(3,3,N_days,Nmc))

  # build matrices
  for(m in seq_along(ua$phi.a.mu)){
  
    phi.11 <- inv.logit(ua$phi.l.mu[m] + ua$beta.11[m]*met)
    phi.22 <- inv.logit(ua$phi.n.mu[m] + ua$beta.22[m]*met)
    lambda <- ifelse(gdd >= 1500 & gdd <= 2500, ua$repro.mu[m], 0)
    theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, inv.logit(ua$grow.na.mu[m]), 0)
    theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)  
    
    # draw transition matrix
    A[1,1,,m] <- phi.11*(1-theta.21)
    A[2,1,,m] <- phi.11*theta.21
    A[2,2,,m] <- phi.22*(1-theta.32)
    A[3,2,,m] <- phi.22*theta.32
    A[3,3,,m] <- inv.logit(ua$phi.a.mu[m] + ua$beta.33[m]*met)
    A[1,3,,m] <- lambda
  }
  return(A)
}