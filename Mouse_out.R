library(ecoforecastR)
library(boot)

load("../FinalOut/GreenControlMR/GreenControlMR_1_5.RData")

source("Functions/model_performance_metrics.R")
source("Functions/mouse_data_jags.R")
data <- suppressWarnings(mouse_data_jags("Green Control"))

days <- data$days
dt <- data$dt
ind <- data$ind
time <- data$time
precip <- data$precip
temp <- data$temp
temp[data$temp.mis] <- mean(temp, na.rm = TRUE)

params <- as.matrix(out$params)
jags.out <- as.matrix(jags.out)
tot.mice.pred <- jags.out[,grep("N[", colnames(jags.out), fixed = TRUE)]
gamma <- jags.out[,grep("gamma[", colnames(jags.out), fixed = TRUE)]

Nmc <- 40
draw <- sample.int(nrow(params), Nmc, replace = TRUE)

lambda.mean <- params[draw, "lambda.mean"]
beta.1 <- params[draw, "beta[1]"]
beta.2 <- params[draw, "beta[2]"]
theta <- params[draw, "theta"]
gamma <- gamma[draw,]

state <- tot.mice.pred[draw,]

X <- suppressWarnings(known_states(data$y))
#X <- X[1:(nrow(X)*0.75),]
Q <- 1-X

# y <- p <- array(0, dim = c(nrow(X), ncol(X), Nmc))
y <- p <- array(0, dim = c(ind, days, Nmc))

# build latent state capture matrix for every day (backfilling from known states matrix)
x.latent <- q <- matrix(0, ind, days)
for(i in seq_along(dt)){
  col.index <- dt[i]
  x.latent[,col.index] <- X[,i]
}
x.latent <- suppressWarnings(known_states(x.latent))
q.latent <- 1 - x.latent

# latent.state <- state.obs <- matrix(0, Nmc, time)
latent.state <- state.obs <- matrix(0, Nmc, days)

for(m in 1:Nmc){
  if(m%%10==0){cat(round(m/Nmc*100), "% ensembles completed\n")}
  
  lambda <- inv.logit(lambda.mean[m] + beta.1[m]*precip[,1] + beta.2[m]*temp[,1])
  
  # phi <- rep(NA, time-1)
  # for(d in 1:(time-1)){  ## loop over capture event
  #   phi[d] <- prod(lambda[dt[d]:(dt[d+1])])
  # } # t
  
  mu1.save <- matrix(NA, time, ind)
  for(d in 2:days){
    ## State Process
    
    t <- max(which(d >= dt)) # code to go from dt to days (for gamma)
    
    q[,d-1] <- 1 - p[,d-1,m]
    # mu1 <- phi[t-1] * X[, t-1] + gamma[m,t] * prod(Q[, 1:(t-1)])
    mu1 <- lambda[d-1] * p[,d-1,m] + gamma[m,t] * prod(q[,1:(d-1)])
    
    # ic.sim is a vector of states (1 or 0) 
    ic.sim <- rbinom(ind, 1, mu1)
    p[, d, m] <- ic.sim
    
    # find which ensembles match the known state from d-1 to d
    
    ## Observation process
    mu2 <- theta[m] * ic.sim
    y[, d, m] <- rbinom(ind, 1, mu2)
    
    
    if(d%%100==0){cat(d, "days completed in ensemble",m,"\n")}
  } # t

  
  latent.state[m,] <- apply(p[,,m], 2, sum) 
  state.obs[m,] <- apply(y[,,m], 2, sum) 
  
}
cat("Finished MCMC Simulation\n")

save(latent.state, state.obs, file = "ExploreScripts/Mice_Pred.RData")
cat("Latent Population Size Saved\n")

latent.state.ens <- p
state.obs.ens <- y
save(latent.state.ens, state.obs.ens, file = "ExploreScripts/Mice_Ens_Pred.RData")
cat("Ensembles Saved\n")


load("ExploreScripts/Mice_Pred.RData")
dim(latent.state)

ci.latent <- apply(latent.state, 2, quantile, c(0.025, 0.5, 0.975))
ci.obs <- apply(state.obs, 2, quantile, c(0.025, 0.5, 0.975))
mna <- apply(X, 2, sum)

time.plot <- 1:ncol(ci.latent)
plot(time.plot, ci.latent[2,], pch = "", ylim = c(0,150))
ciEnvelope(time.plot, ci.latent[1,], ci.latent[3,], col = "lightblue")
lines(time.plot, ci.latent[2,])
points(time.plot, mna, pch = 16)

plot(time.plot, ci.obs[2,], pch = "", ylim = c(0,150))
ciEnvelope(time.plot, ci.obs[1,], ci.obs[3,], col = "lightblue")
lines(time.plot, ci.obs[2,])
points(time.plot, mna, pch = 16)

plot(ci.latent[2,], mna)
abline(0,1)
plot(ci.obs[2,], mna)
abline(0,1)

r_square_one2one(ci.latent[2,], mna)
r_square_one2one(ci.obs[2,], mna)

# mna <- apply(X, 2, sum)
# pre <- apply(y, 2, sum)
# lat <- apply(p, 2, sum)
# par(mfrow=c(1,1))
# time.plot <- 1:length(pre)
# plot(time.plot, mna, type="l")
# lines(time.plot, pre, col=2)
# lines(time.plot, lat, col=3)
# lines(time.plot, state[m,], col=4)



