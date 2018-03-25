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

#inits <- function(){list(x = cjs.init.z(capt.hist, first.cap))}

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
  q[t] <- 1 - (1-exp(-lambda[t]*dt[t]))
  theta[t] <- mean.theta
} # t

## survival transitions; 1 = alive, 0 = dead
for (i in 1:ind){
  for (t in 1:(time - 1)){
    surv[1, i, t, 1] <- phi[t]
    surv[0, i, t, 1] <- 0
    surv[1, i, t, 0] <- q[t]
    surv[0, i, t, 0] <- 1

    ## capture probabilities; the first index is capture (1 = yes, 0 = no)
    ## the second index is alive (1) or dead (2) 
    cap[1, i, t, 1] <- theta[t]
    cap[1, i, t, 0] <- 1 - theta[t]
    cap[0, i, t, 1] <- 0
    cap[0, i, t, 0] <- 1
  } # t
} # i



for(i in 1:ind){
  x[i, 1] <-  1
  for (t in 2:time){
    x[i, t] ~ dcat(surv[x[i, t-1], i, t-1,]) # survival process
    y[i, t] ~ dcat(cap[x[i, t], i, t-1,]) # capture process
  } #t
} # i


} # end model "

j.model <- jags.model(file = textConnection(model1),
                      data = data,
                      #inits = inits,
                      n.chains = 3)

jags.out <- coda.samples(model = j.model,
                         variable.names = c("phi", "N"),
                         n.iter = 1000)



plot(jags.out)