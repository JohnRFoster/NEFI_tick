library(fosteR)
library(mvtnorm)
library(ecoforecastR)
library(boot)

source("Functions/Future_Met.R")
source("Functions/life_stage_ci.R")

data <- cary_ticks_met_JAGS(c("Green Control", "Henry Control", "Tea Control"))

load("/projectnb/dietzelab/fosterj/FinalOut/HB_Partial_Temp/Out_Temp_Partial_HB.RData")

params <- window(params, start = 3000)
predict <- window(predict, start = 3000)

params.mat <- as.matrix(params)
ic.mat <- as.matrix(predict)
n.iter <- nrow(params.mat)

Nmc <- 800 
draw <- sample.int(nrow(params.mat), Nmc)
N_est <- 40
df <- 21
threshold <- c(1500,500,2500)

source("Functions/forecast_with_molting_temp.R")
source("Functions/predict_state_one_temp.R")





# Parameter + IC + Process Uncertainty
predict.P.IC.PROC <- predict_state_one_temp(type = c("parameter", "ic", "process"),
                                                 site = "Henry Control",
                                                 params = params.mat,
                                                 ic = ic.mat,
                                                 data = data,
                                                 Nmc = Nmc,
                                                 draw = draw)

forecast.P.IC.PROC <- forecast_with_molting_temp(type = c("parameter", "ic", "process"),
                                             site = "Henry Control",
                                             params = params.mat,
                                             ic = ic.mat,
                                             N_est = N_est,
                                             df = df,
                                             threshold = threshold,
                                             Nmc = Nmc,
                                             draw = draw)

CI.P.IC.PROC.pred <- life.stage.ci(predict.P.IC.PROC, type = "predict")
CI.P.IC.PROC.forecast <- life.stage.ci(forecast.P.IC.PROC, type = "forecast")

plot(1:N_est, CI.P$Larvae[2,], type = "l", xaxt = "n", main = "Larvae", ylim = c(0, 1000))
ciEnvelope(1:N_est, CI.P.IC.PROC$Larvae[1,], CI.P.IC.PROC$Larvae[3,], col = "grey")
ciEnvelope(1:N_est, CI.P.IC$Larvae[1,], CI.P.IC$Larvae[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Larvae[1,], CI.P$Larvae[3,], col = "purple")
lines(1:N_est, CI.P$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Nymph[2,], type = "l", xaxt = "n", main = "Nymph", ylim = c(0,150))
ciEnvelope(1:N_est, CI.P.IC.PROC$Nymph[1,], CI.P.IC.PROC$Nymph[3,], col = "grey")
ciEnvelope(1:N_est, CI.P.IC$Nymph[1,], CI.P.IC$Nymph[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Nymph[1,], CI.P$Nymph[3,], col = "purple")
lines(1:N_est, CI.P$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Adult[2,], type = "l", xaxt = "n", main = "Adult", ylim = c(0,60))
ciEnvelope(1:N_est, CI.P.IC.PROC$Adult[1,], CI.P.IC.PROC$Adult[3,], col = "grey")
ciEnvelope(1:N_est, CI.P.IC$Adult[1,], CI.P.IC$Adult[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Adult[1,], CI.P$Adult[3,], col = "purple")
lines(1:N_est, CI.P$Adult[2,])
axis(1, at = at, labels = ym.seq[at])


# Process Uncertainty
pred.PROC <- forecast_with_molting_temp(type = "process",
                                        site = "Henry Control",
                                        params = params.mat,
                                        ic = ic.mat,
                                        N_est = N_est,
                                        df = df,
                                        threshold = threshold,
                                        Nmc = Nmc,
                                        draw = draw)

CI.PROC <- life.stage.ci(pred.PROC)

plot(1:N_est, CI.PROC$Larvae[2,], type = "l", xaxt = "n", main = "Larvae", ylim = c(0, 2000))
ciEnvelope(1:N_est, CI.PROC$Larvae[1,], CI.PROC$Larvae[3,], col = "grey")
lines(1:N_est, CI.PROC$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.PROC$Nymph[2,], type = "l", xaxt = "n", main = "Nymph", ylim = c(0,150))
ciEnvelope(1:N_est, CI.PROC$Nymph[1,], CI.PROC$Nymph[3,], col = "grey")
lines(1:N_est, CI.PROC$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.PROC$Adult[2,], type = "l", xaxt = "n", main = "Adult", ylim = c(0,50))
ciEnvelope(1:N_est, CI.PROC$Adult[1,], CI.PROC$Adult[3,], col = "grey")
lines(1:N_est, CI.PROC$Adult[2,])
axis(1, at = at, labels = ym.seq[at])
