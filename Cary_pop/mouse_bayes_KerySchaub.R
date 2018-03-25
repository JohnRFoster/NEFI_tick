library(rjags)
library(coda)

## simulate data
capt.hist.1 <- cbind(matrix(rbinom(10 * 5, 1, 0.7), ncol = 10, nrow = 10), matrix(rbinom(10 * 5, 1, 0.2), ncol = 10, nrow = 10)) # data
capt.hist.2 <- cbind(matrix(rbinom(10 * 5, 1, 0.2), ncol = 10, nrow = 10), matrix(rbinom(10 * 5, 1, 0.7), ncol = 10, nrow = 10))

capt.hist <- rbind(capt.hist.1, capt.hist.2)

## simulate dt's (days between trapping events)
dt <- c(1, 30, 1, 28, 1, 28, 1, 21, 1, 34, 1, 26, 1, 33, 150, 1, 30, 1, 28, 1)

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



## function and call to return the index of each row that has the first 1 (capture)
get.first <- function(x) min(which(x!=0))
first.cap <- apply(capt.hist, 1, get.first)

## Function to create a matrix of initial values for the latent state
cjs.init.z <- function(ch,f){ 
  for (i in 1:dim(ch)[1]){ 
    if (sum(ch[i,])==1) next 
    n2 <- max(which(ch[i,]==1)) 
    ch[i,f[i]:n2] <- NA 
  } 
  for (i in 1:dim(ch)[1]){ 
      ch[i,1:f[i]] <- NA 
  } 
  return(ch) 
}


data <- list(y = capt.hist,  dt = dt, first.cap = first.cap, ind = dim(capt.hist)[1], 
             time = dim(capt.hist)[2], x = known.states(capt.hist),
             a_theta = 1, b_theta = 1, 
             a_lambda = 1, b_lambda = 1)

inits <- function(){list(x = cjs.init.z(capt.hist, first.cap))}



model1 = "
model {

# priors


mean.theta ~ dbeta(a_theta, b_theta)

for(t in 1:(time-1)){
  lambda[t] ~ dgamma(a_lambda, b_lambda)
}

# survival and capture probabilites as a function of time
for(t in 1:(time-1)){
  phi[t] <- lambda[t]*exp(-lambda[t]*dt[t])
  theta[t] <- mean.theta
} # t

for(i in 1:ind){

  x[i, first.cap[i]] <- 1

  for (t in (first.cap[i]+1):time){

    ## State Process
    x[i, t] ~ dbern(mu1[i, t])
    mu1[i, t] <- phi[t-1] * x[i, t-1]

    ## Observation process
    y[i, t] ~ dbern(mu2[i, t])
    mu2[i, t] <- theta[t-1] * x[i, t] 

  } # t
} # i

# abundance estimation

# for(i in 1:ind){
#   for(t in 2:time){
#     N.ind[i, t-1] <- equals(x[i, t], 1)
#   } # i
# } # t


for(t in 2:time){
  N[t] <- sum(x[1:ind, t])
} # t 

}"

j.model <- jags.model(file = textConnection(model1),
                      data = data,
                      inits = inits,
                      n.chains = 3)

jags.out <- coda.samples(model = j.model,
                         variable.names = c("phi", "theta", "N"),
                         n.iter = 1000)


## convert to matrix
jags.mat <- as.matrix(jags.out)

jags.df <- as.data.frame(jags.mat)


jags.means <- apply(jags.df, 2, mean)
plot(1:19, jags.means[1:19])
plot(20:38, jags.means[20:38])

