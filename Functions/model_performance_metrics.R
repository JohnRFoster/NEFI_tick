## root mean squared error
rmse <- function(pred, obs){
  n <- length(pred)
  rmse <- sqrt((1/n)*sum((pred - obs)^2))
  return(rmse)
}

## R-squared calculation relative to 1:1 line
r_square_one2one <- function(pred, obs){
  obs.mean <- mean(obs)
  sum.sq.error <- sum((pred - obs)^2)
  diff.obs.mu <- sum((obs - obs.mean)^2)
  r.square <- 1 - (sum.sq.error/diff.obs.mu)
  return(r.square)
}

## Reduced chi-square
reduced_chi_square <- function(pred, obs){
  n <- length(pred)
  red.chi.sq <- (1 / (n-1))*sum(((pred - obs)^2)/var(obs))
  return(red.chi.sq)
}

## Bayes P Value
bayes_p_val <- function(pred, obs){
  bayes.p.val <- rep(NA, nrow(pred))
  for(gg in 1:nrow(pred)){
    cdf <- ecdf(pred[gg,])
    bayes.p.val[gg] <- cdf(obs[gg])  
  }
  return(bayes.p.val)
}


# mean squared prediction error
mspe <- function(pred, obs){
  n <- length(obs)
  diff.sq <- (obs - pred)^2
  mspe <- sum(diff.sq) / n
  return(mspe)
}
