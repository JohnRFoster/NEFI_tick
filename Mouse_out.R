library(ecoforecastR)
library(boot)

source("Functions/mouse_data_jags.R")
source("Functions/mouse_one.R")

site <- "Henry Control"
file <- "Mouse_out_Henry.RData"
load("../FinalOut/HenryControlMR/HenryControlMR_1_6.RData")
data <- suppressWarnings(mouse_data_jags(site))
jags.out <- as.matrix(jags.out)
Nmc <- 5000
draw <- sample.int(nrow(jags.out), Nmc, replace = TRUE)

X <- suppressWarnings(known_states(data$y))
mna <- colSums(X)
obs <- colSums(data$y)
# names(mna) <- NULL
# names(obs) <- NULL

cat("Partitioning parameter Uncertainty \n")
p <- mouse_one(site, jags.out[draw,], c("parameter"), data)

cat("Partitioning initial condition Uncertainty \n")
ic <- mouse_one(site, jags.out[draw,], c("ic"), data)

cat("Partitioning ic and parameter Uncertainty \n")
p.ic <- mouse_one(site, jags.out[draw,], c("ic","parameter"), data)

cat("Partitioning obs, ic, and parameter Uncertainty \n")
o.p.ic <- mouse_one(site, jags.out[draw,], c("observation","ic","parameter"), data)

cat("Partitioning process, ic, and parameter Uncertainty \n")
pr.p.ic <- mouse_one(site, jags.out[draw,], c("process","ic","parameter"), data)

cat("Partitioning all Uncertainty \n")
o.pr.p.ic <- mouse_one(site, jags.out[draw,], c("observation","process","ic","parameter"), data)

save(p, ic, p.ic, o.p.ic, pr.p.ic, o.pr.p.ic,
     file = file)

# load("Mouse_out.RData")
# 
# 
# latent.state <- ic$latent.state
# state.obs <- ic$latent.state.obs
# 
# quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
# quant.p <- apply(p$latent.state, 2, quantile, quants)
# quant.ic <- apply(ic$latent.state, 2, quantile, quants)
# quant.p.ic <- apply(p.ic$latent.state, 2, quantile, quants)
# quant.o.p.ic <- apply(o.p.ic$latent.state, 2, quantile, quants)
# quant.pr.p.ic <- apply(pr.p.ic$latent.state, 2, quantile, quants)
# quant.o.pr.p.ic <- apply(o.pr.p.ic$latent.state, 2, quantile, quants)
# obs.quant.p <- apply(p$latent.state.obs, 2, quantile, quants)
# obs.quant.ic <- apply(ic$latent.state.obs, 2, quantile, quants)
# obs.quant.p.ic <- apply(p.ic$latent.state.obs, 2, quantile, quants)
# obs.quant.o.p.ic <- apply(o.p.ic$latent.state.obs, 2, quantile, quants)
# obs.quant.pr.p.ic <- apply(pr.p.ic$latent.state.obs, 2, quantile, quants)
# obs.quant.o.pr.p.ic <- apply(o.pr.p.ic$latent.state.obs, 2, quantile, quants)
# 
# time.plot <- 1:(ncol(quant.o.p.ic))
# 
# cols <- c('#edf8fb','#b3cde3','#8c96c6','#88419d')
# cols <- c('#a6611a','#dfc27d','#80cdc1','#018571')
# ylab <- "Individuals"
# xlab <- "Time"
# 
# # ci.obs <- apply(state.obs, 2, quantile, quants)
# 
# # parameter
# plot(time.plot, quant.p[5,], ylim = c(0,max(quant.o.pr.p.ic)+5), pch="")
# ciEnvelope(time.plot, quant.p[1,], quant.p[5,], col = cols[4])
# # ciEnvelope(time.plot, obs.quant.p[1,],obs.quant.p[5,], col = "lightgrey")
# lines(time.plot, quant.p[3,])
# # lines(time.plot, obs)
# 
# # ic
# plot(time.plot, quant.ic[5,], ylim = c(0,max(quant.o.pr.p.ic)+5), pch="",
#      xlab = xlab, ylab = ylab)
# ciEnvelope(time.plot, quant.ic[1,], quant.ic[5,], col = cols[4])
# lines(time.plot, quant.ic[3,])
# legend("topright",
#        c("IC"),
#        lty = 1,
#        col = cols[4])
# 
# 
# # parameter + ic
# plot(time.plot, quant.ic[5,], ylim = c(0,max(quant.o.pr.p.ic)+5), pch="",
#      xlab = xlab, ylab = ylab)
# ciEnvelope(time.plot, quant.p.ic[1,], quant.p.ic[5,], col = cols[3])
# ciEnvelope(time.plot, quant.ic[1,], quant.ic[5,], col = cols[4])
# lines(time.plot, quant.ic[3,])
# legend("topright",
#        c("IC", "IC+Parameter"),
#        lty = rep(1,2),
#        col = cols[4:3])
# 
# # parameter + ic + process
# plot(time.plot, quant.ic[5,], ylim = c(0,max(quant.o.pr.p.ic)+5), pch="",
#      xlab = xlab, ylab = ylab)
# ciEnvelope(time.plot, quant.pr.p.ic[1,], quant.pr.p.ic[5,], col = cols[2])
# ciEnvelope(time.plot, quant.p.ic[1,], quant.p.ic[5,], col = cols[3])
# ciEnvelope(time.plot, quant.ic[1,], quant.ic[5,], col = cols[4])
# lines(time.plot, quant.ic[3,])
# legend("topright",
#        c("IC", "IC+Parameter","IC+Parameter+Process"),
#        lty = rep(1,3),
#        col = cols[4:2])
# 
# # parameter + ic + process with predicted observed
# plot(time.plot, quant.ic[5,], ylim = c(0,max(quant.o.pr.p.ic)+5), pch="",
#      xlab = xlab, ylab = ylab)
# ciEnvelope(time.plot, quant.pr.p.ic[1,], quant.pr.p.ic[5,], col = cols[2])
# ciEnvelope(time.plot, quant.p.ic[1,], quant.p.ic[5,], col = cols[3])
# ciEnvelope(time.plot, quant.ic[1,], quant.ic[5,], col = cols[4])
# ciEnvelope(time.plot, obs.quant.o.pr.p.ic[1,],obs.quant.o.pr.p.ic[5,], col = "lightgrey")
# lines(time.plot, quant.ic[3,])
# lines(time.plot, obs, lty = "dashed")
# legend("topright",
#        c("IC", "IC+Parameter","IC+Parameter+Process", "Predicted Observed", "Observed"),
#        lty = c(rep(1,4), 2),
#        col = c(cols[4:2], "lightgrey", 1))
# 
# 
