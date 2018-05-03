
library(rjags)


y <- read.csv("GreenCaptHist.csv") # capture histories
y <- apply(y, 2, as.numeric)

x <- read.csv("KnownStatesGreen.csv") # known states for each individual
x <- apply(x, 2, as.numeric)

dt <- read.csv("dtGreen.csv") # days between trapping events, in days

data <- list(y = y,  
             dt = dt, 
             ind = nrow(y), 
             time = ncol(y),
             a_theta = 1, b_theta = 1,
             a_lambda = 1, b_lambda = 1)

inits <- function(){list(x = x,
                         lambda = runif(1, 0, 1),
                         theta = runif(1, 0, 1))} 

model3 = "
model {

# priors
theta ~ dbeta(a_theta, b_theta)        # prior on capture probability
lambda ~ dgamma(a_lambda, b_lambda)    # prior on survival probability

for(t in 1:(time-1)){
  phi[t] <- exp(-lambda*dt[t,1])       # calculate survival probability as a function of time
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
  ## B[t] <- sum(recruit[1:ind, t])     # number ot entries
} # t
}"

j.model <- jags.model(file = textConnection(model3),
                      data = data,
                      inits = inits,
                      n.chains = 1)





xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

for (samp in 1:500) {
  if (samp == 1){
    jags.out <- coda.samples(model = j.model,
                             variable.names = c("lambda", "theta", "N", "b", "psi", "gamma"),
                             n.iter = 1000,
                             thin = 10)
  } else {
    jags.out <- update.jags(j.model, 1000)
  }
  saveRDS(jags.out, file = paste("Cary_Recruit_Out_", xx, ".", samp, ".rds", sep = ""))
}



