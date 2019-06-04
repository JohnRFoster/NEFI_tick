# root mean squared error
rmse <- function(pred, obs){
  n <- length(pred)
  rmse <- sqrt(mean((pred - obs)^2))
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