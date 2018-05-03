
library(rjags)

y <- read.csv("GreenCaptHist.csv") # capture histories
y <- apply(y, 2, as.numeric)

ones <- matrix(1, nrow = nrow(y), ncol = ncol(y)) # Matrix for initial conditions

dt <- read.csv("dtGreen.csv") # vector of time differences, in days, between capture events

data <- list(y = y,  
             dt = dt, 
             ind = nrow(y), 
             time = ncol(y),
             a_theta = 1, b_theta = 1,
             a_lambda = 1, b_lambda = 1)

inits <- function(){list(x = ones,
                         lambda = runif(1, 0, 1),
                         theta = runif(1, 0, 1))} 

model1 = "
model {

# priors
theta ~ dbeta(a_theta, b_theta)
lambda ~ dgamma(a_lambda, b_lambda)

for(t in 1:(time-1)){
  phi[t] <- exp(-lambda*dt[t,1])
} # t


for(i in 1:ind){

  x[i, 1] ~ dbern(0.5)

  for (t in (2):time){

    ## State Process
    x[i, t] ~ dbern(mu1[i, t])
    mu1[i, t] <- phi[t-1] * x[i, t-1]

    ## Observation process
    y[i, t] ~ dbern(mu2[i, t])
    mu2[i, t] <- theta * x[i, t] 
  } # t
} # i

# abundance estimation
for(t in 2:(time-1)){
  N[t] <- sum(x[, t])
} # t

}"

j.model <- jags.model(file = textConnection(model1),
                      data = data,
                      inits = inits,
                      n.chains = 1)

xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

for (samp in 1:500) {
  if (samp == 1){
    jags.out <- coda.samples(model = j.model,
                           variable.names = c("lambda", "theta", "N"),
                           n.iter = 1000)
  } else {
    jags.out <- update.jags(j.model, 1000)
  }
  saveRDS(jags.out, file = paste("Cary_Null_Out_", xx, ".", samp, ".rds", sep = ""))
}




