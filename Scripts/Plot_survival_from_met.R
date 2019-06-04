library(ecoforecastR)
library(runjags)
library(fosteR)
library(mvtnorm)
library(ggplot2)
library(grid)
library(gridExtra)
library(boot)

source("Functions/predict_state_one_null.R")
source("Functions/predict_state_one_temp.R")
source("Functions/gg_data_ci_partition.R")
source("Functions/model_performance_metrics.R")

life.stage.ci <- function(pred, quants = c(0.025, 0.5, 0.975)){
  time <- 1:(length(pred[1,,1])-1)
  ci <- list()
  larv <- pred[1,,]
  ci[[1]] <- apply(larv[time,], 1, quantile, quants)
  
  nymph <- pred[2,,]
  ci[[2]] <- apply(nymph[time,], 1, quantile, quants)
  
  adult <- pred[3,,]
  ci[[3]] <- apply(adult[time,], 1, quantile, quants)
  
  names(ci) <- c("Larvae","Nymph","Adult")
  return(ci)
}

load("/projectnb/dietzelab/fosterj/FinalOut/HB_Partial_Temp/Out_Temp_Partial_HB.RData")
params <- window(params, start = 3000)
predict <- window(predict, start = 3000)
summary(params)[[1]]

params.mat <- as.matrix(params)
ic.mat <- as.matrix(predict)
n.iter <- nrow(params.mat)

Nmc <- 500 
draw <- sample.int(nrow(params.mat), Nmc)

sites <- c("Green Control","Henry Control","Tea Control")
data <- cary_ticks_met_JAGS(sites)
N_days <- data$N_days
temp <- data$met[,1,]
for(i in 1:3){
  xx <- data$temp.mis[,i]
  temp[xx,i] <- mean(temp, na.rm = T)
}

nymph.sum <- larv.sum <- data.frame()

i=2
x <- seq(-30, 24, by = 0.01)
l.sum <- matrix(NA, Nmc, length(x))
n.sum <- matrix(NA, Nmc, length(x))
for(t in 1:length(x)){
  l.sum[,t] <- inv.logit(params.mat[draw,"phi.l.mu"] + params.mat[draw,"beta.11"]*x[t] + params.mat[draw,paste("alpha.11[", i, "]", sep = "")])
  n.sum[,t] <- inv.logit(params.mat[draw,"phi.n.mu"] + params.mat[draw,"beta.22"]*x[t] + params.mat[draw,paste("alpha.22[", i, "]", sep = "")]) 
}

nymph.surv.quant <- apply(n.sum, 2, quantile, c(0.025, 0.5, 0.975))
larva.surv.quant <- apply(l.sum, 2, quantile, c(0.025, 0.5, 0.975))

par(mfrow = c(2,1))
plot(x, larva.surv.quant[2,], type = "l", xlab = "Temperature", ylab = "Survival Probability",
     main = "Larvae")
ciEnvelope(x, larva.surv.quant[1,], larva.surv.quant[3,], col = "lightblue")
lines(x, larva.surv.quant[2,])

plot(x, nymph.surv.quant[2,], type = "l", xlab = "Temperature", ylab = "Survival Probability",
     main = "Nymph")
ciEnvelope(x, nymph.surv.quant[1,], nymph.surv.quant[3,], col = "lightblue")
lines(x, nymph.surv.quant[2,])
  
  
  
  

