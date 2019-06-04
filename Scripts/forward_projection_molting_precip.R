library(fosteR)
library(mvtnorm)
library(ecoforecastR)
library(boot)

source("Functions/Future_Met.R")
source("Functions/life_stage_ci.R")
source("Functions/get_last_day.R")

load("/projectnb/dietzelab/fosterj/FinalOut/HB_Partial_Precip/Out_Precip_Partial_HB.RData")

params <- window(params, start = 3000)
predict <- window(predict, start = 3000)

params.mat <- as.matrix(params)
ic.mat <- as.matrix(predict)
n.iter <- nrow(params.mat)

Nmc <- 500 
draw <- sample.int(nrow(params.mat), Nmc)
N_est <- 40
df <- 21
threshold <- c(1500,500,2500)

s <- 2

last <- get_last_day()
last.day <- as.Date(last[s,"last.day"])

date.seq <- seq.Date(last.day, by = df, length.out = N_est)
ym.seq <- format(date.seq, "%Y-%m")
at <- seq(1, N_est, length.out = 7)

par(mfrow = c(3,1))

plot(1:N_est, pred.D[1,,], type = "l", xaxt = "n", main = "Larvae")
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, pred.D[2,,], type = "l", xaxt = "n", main = "Nymph")
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, pred.D[3,,], type = "l", xaxt = "n", main = "Adult")
axis(1, at = at, labels = ym.seq[at])

source("Functions/forecast_with_molting_precip.R")
source("Functions/predict_state_one_precip.R")
data <- cary_ticks_met_JAGS(c("Green Control", "Henry Control", "Tea Control"))

# Parameter Uncertainty
forecast.P <- forecast_with_molting_precip(type = "parameter",
                                         site = "Henry Control",
                                         params = params.mat,
                                         ic = ic.mat,
                                         N_est = N_est,
                                         df = df,
                                         threshold = threshold,
                                         Nmc = Nmc,
                                         draw = draw)

predict.P <- predict_state_one_precip(type = "parameter",
                                    site = "Henry Control",
                                    params = params.mat,
                                    ic = ic.mat,
                                    data = data,
                                    Nmc = Nmc,
                                    draw = draw)

CI.P.forecast <- life.stage.ci(forecast.P, type = "forecast")
CI.P.predict <- life.stage.ci(predict.P, type = "predict")

time.predict <- 1:length(CI.P.predict$Larvae[2,])
time.forecast <- (length(CI.P.predict$Larvae[2,])+1):
  (length(CI.P.predict$Larvae[2,])+length(CI.P.forecast$Larvae[2,]))
time.all <- 1:time.forecast[length(time.forecast)]

col.param <- "purple"
col.param.ic <- "lightgreen"
col.process <- "grey"

plot_run <- function(time, ci.state, life.stage,...){
  plot(time, time, pch = "", ylab = "Individuals", xlab = "", main = life.stage, ...)
}

plot_run(time.all, CI.P.predict$Larvae, "Larvae", ylim = c(0,2000))
ciEnvelope(time.predict, CI.P.predict$Larvae[1,], CI.P.predict$Larvae[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Larvae[1,], CI.P.forecast$Larvae[3,], col = col.param)
lines(time.predict, CI.P.predict$Larvae[2,])
lines(time.forecast, CI.P.forecast$Larvae[2,])

plot_run(time.all, CI.P.predict$Nymph, "Nymph", ylim = c(0,100))
ciEnvelope(time.predict, CI.P.predict$Nymph[1,], CI.P.predict$Nymph[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Nymph[1,], CI.P.forecast$Nymph[3,], col = col.param)
lines(time.predict, CI.P.predict$Nymph[2,])
lines(time.forecast, CI.P.forecast$Nymph[2,])

plot_run(time.all, CI.P.predict$Adult, "Adult", ylim = c(0,20))
ciEnvelope(time.predict, CI.P.predict$Adult[1,], CI.P.predict$Adult[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Adult[1,], CI.P.forecast$Adult[3,], col = col.param)
lines(time.predict, CI.P.predict$Adult[2,])
lines(time.forecast, CI.P.forecast$Adult[2,])

# Parameter + IC Uncertainty
forecast.P.IC <- forecast_with_molting_precip(type = c("parameter", "ic"),
                                            site = "Henry Control",
                                            params = params.mat,
                                            ic = ic.mat,
                                            N_est = N_est,
                                            df = df,
                                            threshold = threshold,
                                            Nmc = Nmc,
                                            draw = draw)

predict.P.IC <- predict_state_one_precip(type = c("parameter", "ic"),
                                       site = "Henry Control",
                                       params = params.mat,
                                       ic = ic.mat,
                                       data = data,
                                       Nmc = Nmc,
                                       draw = draw)

CI.P.IC.forecast <- life.stage.ci(forecast.P.IC, type = "forecast")
CI.P.IC.predict <- life.stage.ci(predict.P.IC, type = "predict")

plot_run(time.all, CI.P.predict$Larvae, "Larvae", ylim = c(0,2000))
ciEnvelope(time.predict, CI.P.IC.predict$Larvae[1,], CI.P.IC.predict$Larvae[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Larvae[1,], CI.P.predict$Larvae[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.IC.forecast$Larvae[1,], CI.P.IC.forecast$Larvae[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.forecast$Larvae[1,], CI.P.forecast$Larvae[3,], col = col.param)
lines(time.predict, CI.P.predict$Larvae[2,])
lines(time.forecast, CI.P.forecast$Larvae[2,])

plot_run(time.all, CI.P.predict$Nymph, "Nymph", ylim = c(0,100))
ciEnvelope(time.predict, CI.P.IC.predict$Nymph[1,], CI.P.IC.predict$Nymph[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.IC.forecast$Nymph[1,], CI.P.IC.forecast$Nymph[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Nymph[1,], CI.P.predict$Nymph[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Nymph[1,], CI.P.forecast$Nymph[3,], col = col.param)
lines(time.predict, CI.P.predict$Nymph[2,])
lines(time.forecast, CI.P.forecast$Nymph[2,])

plot_run(time.all, CI.P.predict$Adult, "Adult", ylim = c(0,20))
ciEnvelope(time.predict, CI.P.IC.predict$Adult[1,], CI.P.IC.predict$Adult[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.IC.forecast$Adult[1,], CI.P.IC.forecast$Adult[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Adult[1,], CI.P.predict$Adult[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Adult[1,], CI.P.forecast$Adult[3,], col = col.param)
lines(time.predict, CI.P.predict$Adult[2,])
lines(time.forecast, CI.P.forecast$Adult[2,])


# Parameter + IC + Process Uncertainty
forecast.P.IC.PROC <- forecast_with_molting_precip(type = c("parameter", "ic", "process"),
                                                 site = "Henry Control",
                                                 params = params.mat,
                                                 ic = ic.mat,
                                                 N_est = N_est,
                                                 df = df,
                                                 threshold = threshold,
                                                 Nmc = Nmc,
                                                 draw = draw)

predict.P.IC.PROC <- predict_state_one_precip(type = c("parameter", "ic", "process"),
                                            site = "Henry Control",
                                            params = params.mat,
                                            ic = ic.mat,
                                            data = data,
                                            Nmc = Nmc,
                                            draw = draw)

CI.P.IC.PROC.forecast <- life.stage.ci(forecast.P.IC.PROC, type = "forecast")
CI.P.IC.PROC.predict <- life.stage.ci(predict.P.IC.PROC, type = "predict")

plot_run(time.all, CI.P.predict$Larvae, "Larvae", ylim = c(0,2000))
ciEnvelope(time.predict, CI.P.IC.PROC.predict$Larvae[1,], CI.P.IC.PROC.predict$Larvae[3,], col = col.process)
ciEnvelope(time.forecast, CI.P.IC.PROC.forecast$Larvae[1,], CI.P.IC.PROC.forecast$Larvae[3,], col = col.process)
ciEnvelope(time.predict, CI.P.IC.predict$Larvae[1,], CI.P.IC.predict$Larvae[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.IC.forecast$Larvae[1,], CI.P.IC.forecast$Larvae[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Larvae[1,], CI.P.predict$Larvae[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Larvae[1,], CI.P.forecast$Larvae[3,], col = col.param)
lines(time.predict, CI.P.predict$Larvae[2,])
lines(time.forecast, CI.P.forecast$Larvae[2,])
lines(time.forecast, CI.P.IC.PROC.forecast$Larvae[2,], col = 2)
points(time.predict, data$y[1,time.predict+1,2], pch = pch)
#axis(1, at = at, labels = ym.seq[at])

plot_run(time.all, CI.P.predict$Nymph, "Nymph", ylim = c(0,100))
ciEnvelope(time.predict, CI.P.IC.PROC.predict$Nymph[1,], CI.P.IC.PROC.predict$Nymph[3,], col = col.process)
ciEnvelope(time.forecast, CI.P.IC.PROC.forecast$Nymph[1,], CI.P.IC.PROC.forecast$Nymph[3,], col = col.process)
ciEnvelope(time.predict, CI.P.IC.predict$Nymph[1,], CI.P.IC.predict$Nymph[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.IC.forecast$Nymph[1,], CI.P.IC.forecast$Nymph[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Nymph[1,], CI.P.predict$Nymph[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Nymph[1,], CI.P.forecast$Nymph[3,], col = col.param)
lines(time.predict, CI.P.predict$Nymph[2,])
lines(time.forecast, CI.P.forecast$Nymph[2,])
lines(time.forecast, CI.P.IC.PROC.forecast$Nymph[2,], col = 2)
points(time.predict, data$y[2,time.predict+1,2], pch = pch)
#axis(1, at = at, labels = ym.seq[at])

plot_run(time.all, CI.P.predict$Adult, "Adult", ylim = c(0,50))
ciEnvelope(time.predict, CI.P.IC.PROC.predict$Adult[1,], CI.P.IC.PROC.predict$Adult[3,], col = col.process)
ciEnvelope(time.forecast, CI.P.IC.PROC.forecast$Adult[1,], CI.P.IC.PROC.forecast$Adult[3,], col = col.process)
ciEnvelope(time.predict, CI.P.IC.predict$Adult[1,], CI.P.IC.predict$Adult[3,], col = col.param.ic)
ciEnvelope(time.forecast, CI.P.IC.forecast$Adult[1,], CI.P.IC.forecast$Adult[3,], col = col.param.ic)
ciEnvelope(time.predict, CI.P.predict$Adult[1,], CI.P.predict$Adult[3,], col = col.param)
ciEnvelope(time.forecast, CI.P.forecast$Adult[1,], CI.P.forecast$Adult[3,], col = col.param)
lines(time.predict, CI.P.predict$Adult[2,])
lines(time.forecast, CI.P.forecast$Adult[2,])
lines(time.forecast, CI.P.IC.PROC.forecast$Adult[2,], col = 2)
points(time.predict, data$y[3,time.predict+1,2], pch = pch)
#axis(1, at = at, labels = ym.seq[at])









pred.D <- forecast_with_molting_precip(type = "deterministic",
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     threshold = threshold,
                                     Nmc = 1,
                                     draw = draw)
s <- 2

last <- get_last_day()
last.day <- as.Date(last[s,"last.day"])

date.seq <- seq.Date(last.day, by = df, length.out = N_est)
ym.seq <- format(date.seq, "%Y-%m")
at <- seq(1, N_est, length.out = 7)

par(mfrow = c(3,1))

plot(1:N_est, pred.D[1,,], type = "l", xaxt = "n", main = "Larvae")
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, pred.D[2,,], type = "l", xaxt = "n", main = "Nymph")
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, pred.D[3,,], type = "l", xaxt = "n", main = "Adult")
axis(1, at = at, labels = ym.seq[at])


# Parameter Uncertainty
pred.P <- forecast_with_molting_precip(type = "parameter",
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     threshold = threshold,
                                     Nmc = Nmc,
                                     draw = draw)

CI.P <- life.stage.ci(pred.P)

plot(1:N_est, CI.P$Larvae[2,], type = "l", xaxt = "n", main = "Larvae")
ciEnvelope(1:N_est, CI.P$Larvae[1,], CI.P$Larvae[3,], col = "purple")
lines(1:N_est, CI.P$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Nymph[2,], type = "l", xaxt = "n", main = "Nymph", ylim = c(0,50))
ciEnvelope(1:N_est, CI.P$Nymph[1,], CI.P$Nymph[3,], col = "purple")
lines(1:N_est, CI.P$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Adult[2,], type = "l", xaxt = "n", main = "Adult", ylim = c(0,30))
ciEnvelope(1:N_est, CI.P$Adult[1,], CI.P$Adult[3,], col = "purple")
lines(1:N_est, CI.P$Adult[2,])
axis(1, at = at, labels = ym.seq[at])


# Parameter + IC Uncertainty
pred.P.IC <- forecast_with_molting_precip(type = c("parameter", "ic"),
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     threshold = threshold,
                                     Nmc = Nmc,
                                     draw = draw)

CI.P.IC <- life.stage.ci(pred.P.IC)

plot(1:N_est, CI.P$Larvae[2,], type = "l", xaxt = "n", main = "Larvae", ylim = c(0, 1000))
ciEnvelope(1:N_est, CI.P.IC$Larvae[1,], CI.P.IC$Larvae[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Larvae[1,], CI.P$Larvae[3,], col = "purple")
lines(1:N_est, CI.P$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Nymph[2,], type = "l", xaxt = "n", main = "Nymph", ylim = c(0,50))
ciEnvelope(1:N_est, CI.P.IC$Nymph[1,], CI.P.IC$Nymph[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Nymph[1,], CI.P$Nymph[3,], col = "purple")
lines(1:N_est, CI.P$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P$Adult[2,], type = "l", xaxt = "n", main = "Adult", ylim = c(0,30))
ciEnvelope(1:N_est, CI.P.IC$Adult[1,], CI.P.IC$Adult[3,], col = "lightgreen")
ciEnvelope(1:N_est, CI.P$Adult[1,], CI.P$Adult[3,], col = "purple")
lines(1:N_est, CI.P$Adult[2,])
axis(1, at = at, labels = ym.seq[at])



# Parameter + IC + Process Uncertainty
pred.P.IC.PROC <- forecast_with_molting_precip(type = c("parameter", "ic", "process"),
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     threshold = threshold,
                                     Nmc = Nmc,
                                     draw = draw)

CI.P.IC.PROC <- life.stage.ci(pred.P.IC.PROC, "forecast")

plot(1:N_est, CI.P.IC.PROC$Larvae[2,], type = "l", xaxt = "n", xlab = "", main = "Larvae", ylim = c(0, 2000))
ciEnvelope(1:N_est, CI.P.IC.PROC$Larvae[1,], CI.P.IC.PROC$Larvae[3,], col = "grey")
# ciEnvelope(1:N_est, CI.P.IC$Larvae[1,], CI.P.IC$Larvae[3,], col = "lightgreen")
# ciEnvelope(1:N_est, CI.P$Larvae[1,], CI.P$Larvae[3,], col = "purple")
lines(1:N_est, CI.P.IC.PROC$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P.IC.PROC$Nymph[2,], type = "l", xaxt = "n", xlab = "", main = "Nymph", ylim = c(0,100))
ciEnvelope(1:N_est, CI.P.IC.PROC$Nymph[1,], CI.P.IC.PROC$Nymph[3,], col = "grey")
# ciEnvelope(1:N_est, CI.P.IC$Nymph[1,], CI.P.IC$Nymph[3,], col = "lightgreen")
# ciEnvelope(1:N_est, CI.P$Nymph[1,], CI.P$Nymph[3,], col = "purple")
lines(1:N_est, CI.P.IC.PROC$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.P.IC.PROC$Adult[2,], type = "l", xaxt = "n", xlab = "", main = "Adult", ylim = c(0,50))
ciEnvelope(1:N_est, CI.P.IC.PROC$Adult[1,], CI.P.IC.PROC$Adult[3,], col = "grey")
# ciEnvelope(1:N_est, CI.P.IC$Adult[1,], CI.P.IC$Adult[3,], col = "lightgreen")
# ciEnvelope(1:N_est, CI.P$Adult[1,], CI.P$Adult[3,], col = "purple")
lines(1:N_est, CI.P.IC.PROC$Adult[2,])
axis(1, at = at, labels = ym.seq[at])



# Process Uncertainty
source("Functions/forecast_with_molting_precip.R")
pred.PROC <- forecast_with_molting_precip(type = "process",
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     threshold = threshold,
                                     Nmc = Nmc,
                                     draw = draw)

CI.PROC <- life.stage.ci(pred.PROC, "forecast")

plot(1:N_est, CI.PROC$Larvae[2,], type = "l", xaxt = "n", xlab = "", main = "Larvae", ylim = c(0, 2000))
ciEnvelope(1:N_est, CI.PROC$Larvae[1,], CI.PROC$Larvae[3,], col = "grey")
lines(1:N_est, CI.PROC$Larvae[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.PROC$Nymph[2,], type = "l", xaxt = "n",  xlab = "", main = "Nymph", ylim = c(0,150))
ciEnvelope(1:N_est, CI.PROC$Nymph[1,], CI.PROC$Nymph[3,], col = "grey")
lines(1:N_est, CI.PROC$Nymph[2,])
axis(1, at = at, labels = ym.seq[at])

plot(1:N_est, CI.PROC$Adult[2,], type = "l", xaxt = "n",  xlab = "", main = "Adult", ylim = c(0,50))
ciEnvelope(1:N_est, CI.PROC$Adult[1,], CI.PROC$Adult[3,], col = "grey")
lines(1:N_est, CI.PROC$Adult[2,])
axis(1, at = at, labels = ym.seq[at])
