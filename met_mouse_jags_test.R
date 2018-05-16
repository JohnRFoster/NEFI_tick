library(rjags)
library(tidyverse)

dt.1 <- c(30, 1, 30, 1, 28, 1, 28, 1, 21, 1, 34, 1, 26, 1, 33, 150, 1, 30, 1, 28, 1)
# write.csv(dt.1, "SimulatedDT.csv")
TIME <- c(1,1+cumsum(dt.1))

# parameters
time <- length(TIME)                  # number of capture occasions
marked <- rep(12, time - 1)           # number newly marke individuals
lambda <- rep(0.0025, time - 1)       # mortality prob
theta <- rep(0.7, time - 1)           # capture prob

## simple version, no recruitment
lambda = lambda[1]
theta = theta[1]
ni = 100
ch <- indiv <- matrix(NA,ni,time)
indiv[,1] <- 1
ch[,1] <- rbinom(ni,1,theta)
for(t in 2:time){
  indiv[,t] <- rbinom(ni,1,exp(-lambda*dt.1[t-1]))*indiv[,t-1]
  ch[,t]  <- rbinom(ni,1,theta)*indiv[,t]
}

known.states <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    if(sum(ch[i,]!=0)){
      n1 <- min(which(ch[i,] == 1))
      n2 <- max(which(ch[i,] == 1))
      state[i, n1:n2] <- 1
    }
  }
  return(state)
}

x <- known.states(ch)
ch <- rbind(ch, matrix(0, ncol=ncol(ch), nrow = nrow(ch)*2))
x <- rbind(x, matrix(0, ncol=ncol(x), nrow = nrow(x)*2))


met <- read.csv("Cary_pop/Met_Cary")

index <- 1:(nrow(met)-sum(dt.1))            # number of possible starting positions for met
first.index <- sample(index, 1)*1           # randomly select the first index
last.index <- first.index + sum(dt.1) - 1   # last index

# precip sequence
precip <- met %>% 
  select(TOT_PREC) %>% 
  slice(first.index:last.index)

#temp sequence 
temp.max <- met %>% 
  select(MAX_TEMP) %>% 
  slice(first.index:last.index)


data <- list(y = ch,  
             dt = dt.1, 
             ind = nrow(ch), 
             time = ncol(ch),
             precip = precip,
             temp = temp.max,
             days = nrow(precip),
             a_theta = 1, b_theta = 1,
             a_lambda = 1, b_lambda = 1)

inits <- function(){list(x = x,
                         lambda = runif(1, 0, 1),
                         theta = runif(1, 0, 1))} 

model3 = "
model {

# priors
theta ~ dbeta(a_theta, b_theta)             # prior on capture probability
lambda.mean ~ dgamma(a_lambda, b_lambda)    # prior on mean hazard rate

for(t in 1:days){
  lambda[t] <- lambda.mean
}

for(t in 1:(time-1)){
  phi[t] <- exp(-lambda*dt[t])         # calculate survival probability as a function of time
} # t                                    

for(t in 1:time){
  gamma[t] ~ dunif(0, 1)               # prior on removal entry probability
}

for(i in 1:ind){

  ## first state process
  x[i, 1] ~ dbern(gamma[1])  # prior on time one for each individual, based on removal entry prob
  mu[i] <- x[i, 1] * theta   # mu = state * capture prob

  ## first observation
  y[i, 1] ~ dbern(mu[i])     # first observation a bernouli trial with prob. mu

  for (t in 2:time){

  ## State Process
  q[i, t-1] <- 1 - x[i, t-1]
  x[i, t] ~ dbern(mu1[i, t])
  mu1[i, t] <- phi[t-1] * x[i, t-1] + gamma[t] * prod(q[i, 1:(t-1)])

  ## Observation process
  y[i, t] ~ dbern(mu2[i, t])
  mu2[i, t] <- theta * x[i, t] 
  } # t
} # i

# calculate derived population parameters
for(t in 1:time){
  qgamma[t] <- 1 - gamma[t]
} # t
cprob[1] <- gamma[1]

for(t in 2:time){
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
} # t

psi <- sum(cprob[])               # inclusion probability

for(t in 1:time){
  b[t] <- cprob[t] / psi          # entry probability
} # t

for(i in 1:ind){
  recruit[i, 1] <- x[i, 1]
  
  for(t in 2:time){
  recruit[i, t] <- (1 - x[i, t-1]) * x[i, t]
  } # t
} # i

for (t in 1:time){
  N[t] <- sum(x[1:ind, t])             # actual population size
  B[t] <- sum(recruit[1:ind, t])     # number ot entries
} # t
}"

j.model <- jags.model(file = textConnection(model3),
                      data = data,
                      inits = inits,
                      n.chains = 3)

jags.out <- coda.samples(model = j.model,
                           variable.names = c("lambda", "theta", "N", "b", "psi", "gamma"),
                           n.iter = 1000)

