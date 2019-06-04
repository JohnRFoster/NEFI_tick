library(ecoforecastR)
library(fosteR)
library(mvtnorm)

sites <- c("Green Control","Henry Control","Tea Control")
data <- cary_ticks_met_JAGS(sites)
data$mna <- mna_tick_current()

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
    data$mna <- data$mna[1,1:data$N_days]
  } else if(site.run == "Henry Control"){
    data$y <- data$y[,,2]
    data$N_est <- data$N_est[2]
    data$N_days <- data$N_days[2]
    data$dt.index <- data$dt.index[2,]
    data$df <- data$df[2,]
    data$met <- data$met[1:data$N_days,,2]
    data$temp.mis <- data$temp.mis[,2]
    data$rh.mis <- data$rh.mis[,2]
    data$mna <- data$mna[2,1:data$N_days]
  } else {
    data$y <- data$y[,-73,3]
    data$N_est <- data$N_est[3]
    data$N_days <- data$N_days[3]
    data$dt.index <- data$dt.index[3,]
    data$df <- data$df[3,]
    data$met <- data$met[1:data$N_days,,3]
    data$temp.mis <- data$temp.mis[,3]
    data$rh.mis <- data$rh.mis[,3]
    data$mna <- data$mna[3,1:data$N_days]
  }
  
  if("temp" %in% met){
    data$temp <- data$met[,1]
  }
  
  if("rh" %in% met){
    data$rh <- data$met[,2]
  }
  
  if("precip" %in% met){
    data$precip <- data$met[,3]
  }
  
  return(data)
}

# load("FinalOut/Independent_Fits/ZIP_Cary_Tick_MNACurrent/GreenControl_MNACurrent_ZIP.RData")
# site.run <- "Green Control"
# s <- 1

load("FinalOut/Independent_Fits/ZIP_Cary_Tick_MNACurrent/HenryControl_MNACurrent_ZIP.RData")
site.run <- "Henry Control"
s <- 2
 
# load("FinalOut/Independent_Fits/ZIP_Cary_Tick_MNACurrent/TeaControl_MNACurrent_ZIP.RData")
# site.run <- "Tea Control"
# s <- 3

data <- data.site.run(data = data, site.run = site.run)
mna <- data$mna
params <- as.matrix(out$params)
ic <- as.matrix(out$predict)
Nmc <- 500
draw <- sample.int(nrow(params), Nmc)

l.phi <- n.phi <- a.phi <- repro <- ln.psi <- na.psi <- matrix(NA,Nmc,length(mna))

for(g in 1:Nmc){
  ln.psi[g,] <- logit(params[draw[g], "grow.ln.mu"] + params[draw[g], "beta.21"]*mna)
  na.psi[g,] <- logit(params[draw[g], "grow.na.mu"] + params[draw[g], "beta.32"]*mna)
}

# these demographic rates are constant, so each is a vector of randomly sampled iterations
repro <- exp(params[draw, "repro.mu"])
l.phi <- logit(params[draw, "phi.l.mu"])
n.phi <- logit(params[draw, "phi.n.mu"])
a.phi <- logit(params[draw, "phi.a.mu"])

l.surv.day <- quantile(l.phi, c(0, 0.025, 0.5, 0.975, 1))
n.surv.day <- quantile(n.phi, c(0, 0.025, 0.5, 0.975, 1))
a.surv.day <- quantile(a.phi, c(0, 0.025, 0.5, 0.975, 1))
repro.day <- quantile(repro, c(0, 0.025, 0.5, 0.975, 1))
ln.day <- apply(ln.psi, 2, quantile, c(0, 0.025, 0.5, 0.975, 1))
na.day <- apply(na.psi, 2, quantile, c(0, 0.025, 0.5, 0.975, 1))

larv.check <- nymph.check <- matrix(NA,Nmc,length(mna))
for(t in 1:data$N_days){ # loop over all days
  larv.check[,t] <- l.phi + ln.psi[,t] 
  nymph.check[,t] <- n.phi + na.psi[,t] 
}

length(which(larv.check > 1)) / (Nmc*length(mna))
length(which(nymph.check > 1)) / (Nmc*length(mna))
  
time <- 1:data.site.run$N_days
# plot(time, larv.check[1,])
# for(g in 2:Nmc){
#   points(time, larv.check[g,])
# }
