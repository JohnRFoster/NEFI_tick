source("Functions/ua_parts.R")
source("Functions/build_periodic_matrices_null.R")
source("Functions/build_periodic_matrices_1proc.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/site_data_met.R")
source("Functions/hb_restructure.R")
source("Functions/sensitivity_elasticity.R")

library(ecoforecastR)

load("../FinalOut/A_Correct/ObsProcModels/Obs_L1N0A0.Proc_AdultRh/Henry/Combined_thinMat_Obs_L1.N0.A0_Proc_ARh_HenryControl.RData")
site <- "Henry Control"
met.variable <- "vpd"

Nmc <- 100
draw <- sample.int(nrow(params.mat), Nmc)
params.mat <- params.mat[draw,]
predict.mat <- predict.mat[draw,]

data.full <- cary_ticks_met_JAGS()
data.site <- site_data_met(site, met.variable, data.full)
N_days <- data.site$N_days
gdd <- data.site$gdd[1:N_days]

ua <- ua_parts(params.mat, predict.mat, c("ic", "parameter", "process"))
A <- build_periodic_matrices_1proc(ua, data.site)

sens <- array(NA, dim = c(6, N_days, Nmc))
sens.x <- matrix(NA, Nmc, N_days)
delta <- 0.01
for(m in 1:Nmc){
  
  phi.11 <- inv.logit(ua$phi.l.mu[m])
  phi.22 <- inv.logit(ua$phi.n.mu[m])
  phi.33 <- inv.logit(ua$phi.a.mu[m])
  theta.32 <- ifelse(gdd <= 1000 | gdd >= 2500, inv.logit(ua$grow.na.mu[m]), 0)
  theta.21 <- ifelse(gdd >= 500 & gdd <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)
  repro <- ua$repro.mu[m]
  
  for(t in 1:N_days){
    # right eigenvector
    w <- eigen(A[,,t,m])$vectors[,1] 
    
    # left eigenvector
    v <- eigen(t(A[,,t,m]))$vectors[,1] 
    
    # not used in calculation, but order of vector matters for building jacobian
    theta.vec <- c(phi.11, phi.22, theta.21[t], theta.32[t], repro, phi.33)
    
    # build jacobian matrix
    # derivative of vec(A[i]) wrt theta.vec[j]
    dAdtheta <- matrix(0, length(as.vector(A[,,t,m])), length(theta.vec))
    dAdtheta[1,1] <- 1-theta.21[t]
    dAdtheta[4,1] <- theta.21[t]
    dAdtheta[5,2] <- 1-theta.32[t]
    dAdtheta[6,2] <- theta.32[t]
    dAdtheta[1,3] <- -1*phi.11  
    dAdtheta[4,3] <- phi.11  
    dAdtheta[5,4] <- -phi.22  
    dAdtheta[6,4] <- -1*phi.22  
    dAdtheta[7,5] <- 1  
    dAdtheta[9,6] <- 1  
    
    # kronecker product
    numerator <- kronecker(w, v)
    
    # row vector v %*% column vector w
    denominator <- v %*% t(t(w))
    
    # sensetivity to lower level parameters for every A (every day)
    sens[,t,m] <- (numerator/as.vector(denominator)) %*% dAdtheta
    sens.x[m,t] <- (numerator/as.vector(denominator)) %*% dAdtheta %*% c(rep(0,5),ua$beta.a[m])
  }
}

# plot(time, nymph.surv.sens.quant[5,], pch="")
# ciEnvelope(time, nymph.surv.sens.quant[1,], nymph.surv.sens.quant[5,], col = "lightblue")
# lines(time, nymph.surv.sens.quant[3,])

quants <- c(0.025,0.25,0.5,0.75,0.975)
time <- 1:N_days
param.names <- c("phi.11", "phi.22", "theta.21", "theta.32", "repro", "phi.33")

jpeg(filename = "Plots/Sensitivity_VPD_LowerLevel.jpeg", width = 12, height = 12)
par(mfrow = c(2,3))
for(j in 1:6){
  nymp.calc <- apply(sens[j,,], 1, quantile, quants)
  # gg.data <- as.data.frame(t(nymp.calc))
  # names(gg.data) <- c("low", "q1", "median", "q3", "high")
  # gg.data$time <- time
  # 
  # p <- ggplot(gg.data) +
  #        geom_ribbon(aes(x = time, ymin = low, ymax = high), color = "lightblue", alpha = 0.3) +
  #       # geom_line(aes(time, median)) +
  #        theme_classic()
  # 
  # print(p)
  
  plot(time, nymp.calc[5,], pch="", main = param.names[j], ylim = c(range(nymp.calc)))
  ciEnvelope(time, nymp.calc[1,], nymp.calc[5,], col = "lightblue")
  lines(time, nymp.calc[3,])  
}
# mtext("Sensitivity to Lower Leverl Params: VPD", outer = TRUE, cex = 1.5)
dev.off()

par(mfrow = c(1,1))
x.quantile <- apply(sens.x, 2, quantile, quants)
plot(time, x.quantile[5,], pch="", main = "Sensitivity to x", ylim = c(range(x.quantile)))
ciEnvelope(time, x.quantile[1,], x.quantile[5,], col = "lightblue")
lines(time, x.quantile[3,])  








# w <- v <- array(NA, dim = c(3,dim(A)[3],Nmc))
# theta.21 <- theta.32 <- matrix(NA, Nmc, N_days)
# for(m in 1:Nmc){
#   for(t in 1:N_days){
#     
#     theta.32[m,t] <- ifelse(gdd[t] <= 1000 | gdd[t] >= 2500, inv.logit(ua$grow.na.mu[m]), 0)
#     theta.21[m,t] <- ifelse(gdd[t] >= 500 & gdd[t] <= 2500, inv.logit(ua$grow.ln.mu[m]), 0)
#   }
# }




# theta.21.sens <- theta.32.sens <- phi.11.sens <- phi.22.sens <- matrix(NA, Nmc, N_days)
# for(m in 1:Nmc){
#   phi.11 <- inv.logit(ua$phi.l.mu[m])
#   phi.22 <- inv.logit(ua$phi.n.mu[m])
#   for(t in 1:N_days){
#     denominator <- v[,t,m] %*% t(t(w[,t,m]))
#     # larvae
#     phi.11.sens[m,t] <- (w[1,t,m]*(v[1,t,m]+(theta.21[m,t]*(v[2,t,m]-v[1,t,m])))) / denominator
#     theta.21.sens[m,t] <- (phi.11*w[1,t,m]*(v[2,t,m]-v[1,t,m])) / denominator
#     
#     # nymph
#     phi.22.sens[m,t] <- (w[2,t,m]*(v[2,t,m]+(theta.32[m,t]*(v[3,t,m]-v[2,t,m])))) / denominator
#     theta.32.sens[m,t] <- (phi.22*w[1,t,m]*(v[2,t,m]-v[1,t,m])) / denominator
#     
#   }
# }
# 
# 
# larva.surv.sens.quant <- apply(phi.11.sens, 2, quantile, quants)
# larva.tran.sens.quant <- apply(theta.21.sens, 2, quantile, quants)
# nymph.surv.sens.quant <- apply(phi.22.sens, 2, quantile, quants)
# nymph.tran.sens.quant <- apply(theta.32.sens, 2, quantile, quants)
# time <- 1:ncol(larva.surv.sens.quant)