library(ecoforecastR)

source("Functions/ammend_chains.R")
source("Functions/split_mice.R")
source("Functions/plot_posterior_chains.R")
source("Functions/plot_post_prior.R")

dir <- "../FinalOut/HB_Partial_GDD/Mice/WindowLoop/GDDSwitch_Mice_Win"
num.chains <- 5
num.out <- 8

out <- ammend_chains(dir, num.chains, num.out, save = "GDDSwitch_Mice_Win")

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



dir <- "../FinalOut/Independent_Fits/GDDThreshold/Temp_Obs_Correct/K_set/LarvaOnly/beta_1/"
# site.folders <- c("Green", "Henry", "Tea")
site.folders <- c("Henry", "Tea")
model <- c(#"Combined_thinMat_Temp_obs_beta_2_K_set_GreenControl.RData",
           "Combined_thinMat_Obs_LarvaOnly_Beta1_HenryControl.RData",
           "Combined_thinMat_Obs_LarvaOnly_Beta1_TeaControl.RData")

beta.1 <- list()
for(i in 1:length(model)){
  load(file.path(dir, site.folders[i], model[i]))
  # out$params <- window(out$params, 2500)
  beta.1[[i]] <- params.mat
  # print(model[i])
  # print(summary(out$params))
  # rm(out)
}

source("Functions/plot_site_posteriors.R")

all.params <- rbind(as.matrix(beta.1[[1]]),
                    as.matrix(beta.1[[2]]))
                    # as.matrix(beta.1[[3]]))

plot_site_posteriors(all.params[,-c(1:9)], n.site = 2)
