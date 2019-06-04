library(ecoforecastR)
library(fosteR)
library(mvtnorm)

sites <- c("Green Control","Henry Control","Tea Control")
data <- cary_ticks_met_JAGS(sites)

# load("FinalOut/Independent_Fits/ZIP_Cary_Tick_RH/GreenControl_RH_ZIP.RData")
# site.run <- "Green Control"
# s <- 1

# not converged
# load("FinalOut/Independent_Fits/Poisson_Cary_Tick_RH/HenryControl_RH_Pois.RData")

load("FinalOut/Independent_Fits/ZIP_Cary_Tick_RH/HenryControl_RH_ZIP.RData")
site.run <- "Henry Control"
s <- 2
 
# load("FinalOut/Independent_Fits/ZIP_Cary_Tick_RH/TeaControl_RH_ZIP.RData")
# site.run <- "Tea Control"
# s <- 3

params <- as.matrix(out$params)
params <- params[1000:5000,]
ic <- as.matrix(out$predict)
ic <- ic[1000:5000,]
Nmc <- 500
draw <- sample.int(nrow(params), Nmc)

data.site.run <- function(data, site.run, met = NULL){
  if(site.run == "Green Control"){
    data$y <- data$y[,-73,1]
    data$N_est <- data$N_est[1]
    data$N_days <- data$N_days[1]
    data$dt.index <- data$dt.index[1,]
    data$df <- data$df[1,]
    data$met <- data$met[1:data$N_days,,1]
    data$temp.mis <- data$temp.mis[,1]
    data$rh.mis <- data$rh.mis[,1]
  } else if(site.run == "Henry Control"){
    data$y <- data$y[,,2]
    data$N_est <- data$N_est[2]
    data$N_days <- data$N_days[2]
    data$dt.index <- data$dt.index[2,]
    data$df <- data$df[2,]
    data$met <- data$met[1:data$N_days,,2]
    data$temp.mis <- data$temp.mis[,2]
    data$rh.mis <- data$rh.mis[,2]
  } else {
    data$y <- data$y[,-73,3]
    data$N_est <- data$N_est[3]
    data$N_days <- data$N_days[3]
    data$dt.index <- data$dt.index[3,]
    data$df <- data$df[3,]
    data$met <- data$met[1:data$N_days,,3]
    data$temp.mis <- data$temp.mis[,3]
    data$rh.mis <- data$rh.mis[,3]
  }
  
  if("temp" %in% met){
    data$temp <- data$met[,1]
    mean.temp <- mean(data$temp, na.rm = TRUE)
    for(i in 1:length(data$temp.mis)){
      data$temp[data$temp.mis[i]] <- mean.temp
    }
  }
  
  if("rh" %in% met){
    data$rh <- data$met[,2]
    mean.rh <- mean(data$rh, na.rm = TRUE)
    for(i in 1:length(data$rh.mis)){
      data$rh[data$rh.mis[i]] <- mean.rh
    }
  }
  
  if("precip" %in% met){
    data$precip <- data$met[,3]
  }
  
  return(data)
}

data.site.run <- data.site.run(data = data, site.run = site.run, met = "rh")

rh <- data.site.run$rh

l.phi <- n.phi <- a.phi <- repro <- ln.psi <- na.psi <- matrix(NA,Nmc,length(rh))

for(g in 1:Nmc){
  l.phi[g,] <- logit(params[draw[g], "phi.l.mu"] + params[draw[g], "beta.11"]*rh)
  n.phi[g,] <- logit(params[draw[g], "phi.n.mu"] + params[draw[g], "beta.22"]*rh)
  a.phi[g,] <- logit(params[draw[g], "phi.a.mu"] + params[draw[g], "beta.33"]*rh)
}

# these demographic rates are constant, so each is a vector of randomly sampled iterations
repro <- exp(params[draw, "repro.mu"])
ln.psi <- logit(params[draw, "grow.ln.mu"])
na.psi <- logit(params[draw, "grow.na.mu"])

l.surv.day <- apply(l.phi, 2, quantile, c(0, 0.025 ,0.5, 0.95, 1))
n.surv.day <- apply(n.phi, 2, quantile, c(0, 0.025, 0.5, 0.95, 1))
a.surv.day <- apply(a.phi, 2, quantile, c(0, 0.025, 0.5, 0.95, 1))
repro.day <- quantile(repro, c(0, 0.025, 0.5, 0.95, 1))
ln.day <- quantile(ln.psi, c(0, 0.025, 0.5, 0.95, 1))
na.day <- quantile(na.psi, c(0, 0.025, 0.5, 0.95, 1))

larv.check <- nymph.check <- matrix(NA,Nmc,length(rh))
for(t in 1:data.site.run$N_days){ # loop over all days
  larv.check[,t] <- l.phi[,t] + ln.psi 
  nymph.check[,t] <- n.phi[,t] + na.psi 
}

length(which(larv.check > 1)) / (Nmc*length(rh))
length(which(nymph.check > 1)) / (Nmc*length(rh))
  
time <- 1:data.site.run$N_days
# plot(time, larv.check[1,])
# for(g in 2:Nmc){
#   points(time, larv.check[g,])
# }
