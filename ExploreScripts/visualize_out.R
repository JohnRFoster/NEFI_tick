library(ecoforecastR)

source("Functions/ammend_chains.R")
source("Functions/split_mice.R")
source("Functions/plot_posterior_chains.R")
source("Functions/plot_post_prior.R")

dir <- "../FinalOut/HB_Partial_GDD/K_estimate/WindowLoop_4alpha/GDDSwitch_K_Window4alpha"
num.chains <- 5
num.out <- 19

out <- ammend_chains(dir, num.chains, num.out, save = "GDDSwitch_K_Window4alpha")

plot(out$params)
plot_posterior_chains(out$params)
plot_posterior_chains(out$params, c("alpha.11[3]",
                                    "alpha.13[3]",
                                    "alpha.k0[3]",
                                    "alpha.k1[3]",
                                    "deviance",
                                    "k.l2n.low",
                                    "phi.n.mu",
                                    "repro.mu"))


dir <- "../FinalOut/HB_Partial_GDD/Mice/WindowLoop/GDDSwitch_Mice_Win"
num.chains <- 5
num.out <- 2

out <- ammend_chains(dir, num.chains, num.out)
model <- split_mice(out)

plot(model$params)
plot_posterior_chains(model$params)
plot_posterior_chains(model$params, c("deviance",
                                      "grow.na.mu",
                                      "phi.a.mu",
                                      "phi.l.mu",
                                      "phi.n.mu",
                                      "repro.mu"))

for(i in 1:5){
  print(i)
  print(model$params[[i]][1,"beta.21"])
  print(range(model$params[[i]][,"beta.21"]))
  print(model$params[[i]][1,"beta.32"])
  print(range(model$params[[i]][,"beta.32"]))
}



mice.est <- as.matrix(model$mice)

mice.ci <- apply(mice.est, 2, quantile, c(0.025, 0.5, 0.975))
par(mfrow = c(1,1))

time <- 1:ncol(mice.ci)
plot(time, mice.ci[2,], type = "l")
ciEnvelope(time, mice.ci[1,], mice.ci[3,], col = "lightblue")
lines(time, mice.ci[2,])



dir <- "../FinalOut/Independent_Fits/GDDThreshold/WithMice/WindowLoop"
model <- c("GreenControl_GDDSwitch_Low_Mice_AllChains.RData",
           "HenryControl_GDDSwitch_Low_Mice_AllChains.RData",
           "TeaControl_GDDSwitch_Low_Mice_AllChains.RData")

mice <- list()
for(i in 1:3){
  load(file.path(dir, model[i]))
  # out$params <- window(out$params, 2500)
  mice[[i]] <- out$params
  print(model[i])
  # print(summary(out$params))
  rm(out)
}

source("Functions/plot_site_posteriors.R")

all.params <- rbind(as.matrix(mice[[1]]),
                    as.matrix(mice[[2]]),
                    as.matrix(mice[[3]]))

plot_site_posteriors(all.params[,-c(1:9)])
