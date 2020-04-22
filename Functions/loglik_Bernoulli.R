loglik_Bernoulli = function(theta, x) {
  # theta   success probability parameter
  # x       vector of data
  n <- length(x)
  ans <- log(theta)*sum(x)+log(1-theta)*(n-sum(x))
  return(ans)
}