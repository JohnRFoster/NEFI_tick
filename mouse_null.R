library(rjags)
library(coda)

## simulate data ##
# parameters
time <- 20                  # number of capture occasions
marked <- rep(12, time - 1) # number newly marke individuals
phi <- rep(0.65, time - 1)  # survival prob
theta <- rep(0.4, time - 1) # capture prob

# matrix with survival and recapture probabilities
phi.mtx <- matrix(phi, ncol = time - 1, nrow = sum(marked))
theta.mtx <- matrix(theta, ncol = time - 1, nrow = sum(marked))

# function to simulate capture history matrix
simul.cjs <- function(phi.mtx, theta.mtx, marked){
  time <- dim(phi.mtx)[2] + 1
  ch <- matrix(0, ncol = time, nrow = sum(marked)) 
  # Define a vector with the occasion of marking 
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix 
  for (i in 1:sum(marked)){ 
    ch[i, mark.occ[i]] <- 1 
    # Write an 1 at the release occasion 
    if (mark.occ[i] == time) next 
    for (t in (mark.occ[i] + 1):time){ 
      # Bernoulli trial: does individual survive occasion? 
      sur <- rbinom(1, 1, phi.mtx[i, t - 1]) 
      if (sur == 0) break # If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, theta.mtx[i, t - 1]) 
      if (rp == 1) ch[i, t] <- 1 } #t 
  } #i 
  return(ch) 
}

# capture histories
capt.hist <- simul.cjs(phi.mtx, theta.mtx, marked)
y <- rbind(capt.hist, matrix(0, ncol = ncol(capt.hist), nrow = 150))

## function to create a matrix of known latent states (change 0's to 1's in between 1's, NA otherwise)
known.states <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] == 1))
    n2 <- max(which(ch[i,] == 1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

# known states
x <- known.states(capt.hist)

## function and call to return the index of each row that has the first 1 (capture)
get.first <- function(x) min(which(x!=0))
first.cap <- apply(y, 1, get.first)

## simulate dt's (days between trapping events)
dt <- c(1, 30, 1, 28, 1, 28, 1, 21, 1, 34, 1, 26, 1, 33, 150, 1, 30, 1, 28, 1)


data <- list(y = y,  dt = dt, ind = nrow(y), time = ncol(y), first.cap = first.cap,
             a_theta = 1, b_theta = 1, 
             a_lambda = 1, b_lambda = 1)

inits <- function(){list(x = x)}



model1 = "
model {

# priors
theta ~ dbeta(a_theta, b_theta)

# for(i in 1:ind){
#   x[i, 1] ~ dbern(0.5)
# }

for(t in 1:(time-1)){
  lambda[t] ~ dgamma(a_lambda, b_lambda)
}

for(t in 1:(time-1)){
  phi[t] <- lambda[t]*exp(-lambda[t]*dt[t])
} # t

 
for(i in 1:ind){
  for (t in (first.cap[i]+1):time){

    ## State Process
    x[i, t] ~ dbern(mu1[i, t])
    mu1[i, t] <- phi[t-1] * x[i, t-1]

    ## Observation process
    y[i, t] ~ dbern(mu2[i, t])
    mu2[i, t] <- theta * x[i, t] 

  } # t
} # i


# # abundance estimation
# for(i in 1:ind){
#   for(t in 2:time){
#     Nind[i, t-1] <- equals(x[i, t], 1)
#   } # i
# } # t
# 
for(t in 1:time){
  N[t] <- sum(x[1:ind, t])
}

}"

j.model <- jags.model(file = textConnection(model1),
                      data = data,
                      inits = inits,
                      n.chains = 3)

jags.out <- coda.samples(model = j.model,
                         variable.names = c("lambda", "theta", "N"),
                         n.iter = 100)

plot(jags.out)
## convert to matrix
jags.mat <- as.matrix(jags.out)

jags.df <- as.data.frame(jags.mat)


jags.means <- apply(jags.df, 2, mean)
plot(1:19, jags.means[1:19])


