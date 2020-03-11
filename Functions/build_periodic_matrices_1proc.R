build_periodic_matrices_1proc <- function(params, gdd, N_days){
  
  # pull needed data from list
  N_days <- data$N_days
  gdd <- data$gdd[1:N_days]
  met <- data$met[1:N_days]
  met[data$met.mis] <- mean(met, na.rm = TRUE)
  
  # storage
  A <- array(0, dim=c(3,3,N_days,nrow(params)))

  phi.l.mu <- params[,"phi.l.mu"]
  phi.n.mu <- params[,"phi.n.mu"]
  phi.a.mu <- params[,"phi.a.mu"]
  grow.ln.mu <- params[,"grow.ln.mu"]
  grow.na.mu <- params[,"grow.na.mu"]
  repro.mu <- params[,"repro.mu"]
  beta.11 <- params[, "beta.11"]
  beta.22 <- params[, "beta.22"]
  beta.33 <- params[, "beta.33"]
  
  # build matrices
  for(m in seq_along(phi.a.mu)){
  
    phi.11 <- inv.logit(phi.l.mu[m] + beta.11[m]*met)
    phi.22 <- inv.logit(phi.n.mu[m] + beta.22[m]*met)
    lambda <- ifelse(gdd >= 1500 & gdd <= 2500, repro.mu[m], 0)
    theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, inv.logit(grow.na.mu[m]), 0)
    theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, inv.logit(grow.ln.mu[m]), 0)  
    
    # draw transition matrix
    A[1,1,,m] <- phi.11*(1-theta.21)
    A[2,1,,m] <- phi.11*theta.21
    A[2,2,,m] <- phi.22*(1-theta.32)
    A[3,2,,m] <- phi.22*theta.32
    A[3,3,,m] <- inv.logit(phi.a.mu[m] + beta.33[m]*met)
    A[1,3,,m] <- lambda
  }
  return(A)
}