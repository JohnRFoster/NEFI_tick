## attempting to use a questing matrix that is aggregated to the daily scale
## and seasonal transition matrices that are used at each observation
## trying to fit this model in JAGS




library(rjags)
library(fosteR)
library(ecoforecastR)

sites <- c("Green Control")
out.name <- "GreenControl_ZIP_cont"

## file path to output folder
out.folder <- "FinalOut/Independent_Fits/ZIP_Cary_Tick/"

out.path <- paste(out.folder,out.name, sep = "")

data <- cary_ticks_JAGS(sites)
data$R <- diag(1,3,3)
data$R.spring <- diag(1,2,2)
data$R.summer <- diag(1,2,2)
data$R.fall <- diag(1,2,2)
data$R.winter <- diag(1,2,2)
data$y <- data$y[,,1]

inits <- function(){list(#m = data$y,
                         # summer[1,] = data$y[1,],
                         # spring[1,] = data$y[2,],
                         # fall[2,] = data$y[3,],
                         repro.mu = runif(1, 1, 3),
                         phi.l.mu = rnorm(1, 3.501, 1.4993),
                         phi.n.mu = rnorm(1, 5.0631, 1.024),
                         phi.a.mu = rnorm(1, 5, 1),
                         quest.larv2molt.nymph = rnorm(1, -3.501, 1.4993),
                         quest.nymph2molt.adult = rnorm(1, -5.0631, 1.024),
                         molt.nymph.surv = rnorm(1, 5.0639, 0.001))}

monitor <- c("m", 
             "phi.l.mu",
             "phi.n.mu",
             "phi.a.mu",
             "quest.larv2molt.nymph",
             "quest.nymph2molt.adult",
             "repro.mu",
             "SIGMA",
             "SIGMA.spring",
             "SIGMA.summer",
             "SIGMA.fall",
             "SIGMA.winter",
             "egg2questlarv",
             "molt.adult2quest.adult",
             "molt.nymph.surv",
             "molt.nymph2quest.nymph",
             "theta.larvae",
             "theta.nymph",
             "theta.adult",
             "summer",
             "spring",
             "fall",
             "winter",
             "deviance")


model = " model {

### global hyperpriors
phi.l.mu ~ dnorm(3.7815,1.7488)               # larvae survival
phi.n.mu ~ dnorm(3.6241,0.5005)               # nymph survival
phi.a.mu ~ dnorm(5,1)               # adult survival
quest.larv2molt.nymph ~ dnorm(-3.501,0.01)             # larvae -> nymph transition 
quest.nymph2molt.adult ~ dnorm(-5.0631,0.01)             # nymph -> adult transition
repro.mu ~ dnorm(2,1)              # adult -> larvae transition (reproduction)
egg2questlarv ~ dunif(0.8,1)      
quest.adult2dorm.adult ~ dunif(0.8,1)
molt.adult2quest.adult ~ dunif(0.8,1)
molt.nymph.surv ~ dnorm(5.0639,0.50228)       
molt.nymph2quest.nymph ~ dunif(0.8,1)

### precision priors
SIGMA        ~ dwish(R, 4)         # mvn [3 x 3] site process
SIGMA.summer ~ dwish(R.summer, 3)  # mvn [2 x 2] phenology process
SIGMA.spring ~ dwish(R.spring, 3)  # mvn [2 x 2] phenology process
SIGMA.fall   ~ dwish(R.fall, 3)    # mvn [2 x 2] phenology process
SIGMA.winter ~ dwish(R.winter, 3)  # mvn [2 x 2] phenology process

### prior probablity we observe ticks by life stage
theta.larvae ~ dunif(0,1)
theta.nymph ~ dunif(0,1)
theta.adult ~ dunif(0,1)

### prior on first latent process 
spring[1,1] ~ dpois(1000)     # eggs 
summer[1,1] ~ dpois(750)      # questing larvae
fall[1,1]   ~ dpois(500)      # molting nymphs
winter[1,1] ~ dpois(350)      # molting nymphs
spring[2,1] ~ dpois(250)      # questing nymphs
summer[2,1] ~ dpois(100)      # molting adults
fall[2,1]   ~ dpois(50)       # questing adults
winter[2,1] ~ dpois(25)       # dormant adults
x[1,1] <- summer[1,1]
x[2,1] <- spring[2,1]
x[3,1] <- fall[2,1]

### define parameters
for(t in 1:N_days){   # loop over every day in time series

  ## site specific  daily transtion matrices;
  ## parameters are a linear combination of
  ## the global parameter (intercept) and a 
  ## site random effect
  
  logit(A.day[1,1,t]) <- phi.l.mu
  logit(A.day[2,2,t]) <- phi.n.mu
  logit(A.day[3,3,t]) <- phi.a.mu
  A.day[1,2,t] <- 0
  A.day[2,3,t] <- 0
  A.day[3,1,t] <- 0
  A.day[2,1,t] <- 0
  A.day[3,2,t] <- 0
  A.day[1,3,t] <- 0
}

### aggregate daily matrix between sampling events

for(t in 1:(N_est-1)){  ## number of days to estimate latent state

  ## df is the number of days between each sampling occasion
  ## dt.index the number of days elaplsed from first sampling occasion
  
  for(i in 1){  
    TRANS[1:3,1:3,dt.index[1,t]-df[1,t]+1] <- A.day[1:3,1:3,dt.index[1,t]-df[1,t]+1] 
  }                                                          
  for(i in 2:df[1,t]){
    TRANS[1:3,1:3,dt.index[1,t]-df[1,t]+i] <- TRANS[1:3,1:3,dt.index[1,t]-df[1,t]+i-1] %*% 
    A.day[1:3,1:3,dt.index[1,t]-df[1,t]+i]
  }
}


### Process Model
for(t in 1:(N_est-1)){  # number of days for each site

  S[1,1,t] <- egg2questlarv            # egg to questing larva
  S[2,2,t] <- quest.nymph2molt.adult   # questing nymph to molting adult
  T[1,1,t] <- quest.larv2molt.nymph    # questing larva to molting nymph
  T[2,2,t] <- molt.adult2quest.adult   # molting adult to questing adult
  logit(U[1,1,t]) <- molt.nymph.surv          # molting nymph survival
  U[2,2,t] <- quest.adult2dorm.adult   # questing adult to dormant adult transition
  V[1,1,t] <- molt.nymph2quest.nymph   # molting nymph to questing nymph
  V[2,2,t] <- repro.mu                 # dormant adult to egg (reproduction)

  S[1,2,t] <- 0
  S[2,1,t] <- 0
  T[1,2,t] <- 0
  T[2,1,t] <- 0
  U[1,2,t] <- 0
  U[2,1,t] <- 0
  V[1,2,t] <- 0
  V[2,1,t] <- 0

  ## TRANS is the aggrigated transition matrix between observations
  # Ex.quest[1:3,t] <- TRANS[1:3,1:3,dt.index[1,t]] %*% x[1:3,t]
  Ex.quest[1:3,t] <- TRANS[1:3,1:3,dt.index[1,t]] %*% c(summer[1,t],spring[1,t],fall[2,t])
  quest.x[1:3,t+1] ~ dmnorm(Ex.quest[1:3,t], SIGMA)

  ## predict seasonal transitions
  Ex.summer[1:2,t] <- S[1:2,1:2,t] %*% spring[1:2,t]  # c(e[t], x[2,t])
  Ex.fall[1:2,t]   <- T[1:2,1:2,t] %*% summer[1:2,t]  # c(x[1,t], m.a)
  Ex.winter[1:2,t] <- U[1:2,1:2,t] %*% fall[1:2,t]    # c(m.n, x[3,t])
  Ex.spring[1:2,t] <- V[1:2,1:2,t] %*% winter[1:2,t]  # c(m.n, d.a)
  
  ## seasonal abundance 
  spring.x[1:2,t+1] ~ dmnorm(Ex.spring[1:2,t], SIGMA.spring)
  summer.x[1:2,t+1] ~ dmnorm(Ex.summer[1:2,t], SIGMA.summer)
  fall.x[1:2,t+1]   ~ dmnorm(Ex.fall[1:2,t], SIGMA.fall)
  winter[1:2,t+1] ~ dmnorm(Ex.winter[1:2,t], SIGMA.winter)

  ## new questing
  summer[1,t+1] <- quest.x[1,t+1] + summer.x[1,t+1] 
  spring[1,t+1] <- quest.x[2,t+1] + spring.x[1,t+1]
  fall[2,t+1] <- quest.x[3,t+1] + fall.x[2,t+1]

  # summer[1,t+1] <- x[1,t+1]
  summer[2,t+1] <- summer.x[2,t+1]
  # spring[1,t+1] <- x[2,t+1]
  spring[2,t+1] <- spring.x[2,t+1]
  # fall[2,t+1]   <- x[3,t+1]
  fall[1,t+1]   <- fall.x[1,t+1]
  
## Data Model ##
  
  ## fit the blended model to observed data 
  y[1,t] ~ dpois(m[1,t])
  y[2,t] ~ dpois(m[2,t])
  y[3,t] ~ dpois(m[3,t])
  
  ## blend the poisson and zero inflation models
  # m[1,t] <- x[1,t]*b.larvae[t] + 1E-10
  # m[2,t] <- x[2,t]*b.nymph[t] + 1E-10
  # m[3,t] <- x[3,t]*b.adult[t] + 1E-10
  m[1,t] <- summer[1,t]*b.larvae[t] + 1E-10
  m[2,t] <- spring[1,t]*b.nymph[t] + 1E-10
  m[3,t] <- fall[2,t]*b.adult[t] + 1E-10
  
  ## binary outcome of observation by life stage
  b.larvae[t] ~ dbern(theta.larvae)
  b.nymph[t] ~ dbern(theta.nymph)
  b.adult[t] ~ dbern(theta.adult)

} # t

}" # end model

load.module("glm")
load.module("dic")
j.model <- jags.model(file = textConnection(model),
                      data = data,
                      inits = inits,
                    #  n.adapt = 100000,
                      n.chains = 3)
n.iter <- 1000
jags.out <- coda.samples(model = j.model,
                         variable.names = monitor,
                        # burnin = 50000,
                        # thin = 100,
                         n.iter = n.iter)
## split output
out <- list(params = NULL, m = NULL, spring = NULL, summer = NULL, fall = NULL, winter = NULL)
mfit <- as.matrix(jags.out, chains = TRUE)
m.cols <- grep("m[", colnames(mfit), fixed = TRUE)
spring.cols <- grep("spring[", colnames(mfit), fixed = TRUE)
summer.cols <- grep("summer[", colnames(mfit), fixed = TRUE)
fall.cols <- grep("fall[", colnames(mfit), fixed = TRUE)
winter.cols <- grep("winter[", colnames(mfit), fixed = TRUE)
chain.col <- which(colnames(mfit) == "CHAIN")
out$m <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, m.cols)])
out$spring <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, spring.cols)])
out$summer <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, summer.cols)])
out$fall <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, fall.cols)])
out$winter <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, winter.cols)])
out$params <- ecoforecastR::mat2mcmc.list(mfit[, -c(m.cols,spring.cols,summer.cols,fall.cols,winter.cols)])
