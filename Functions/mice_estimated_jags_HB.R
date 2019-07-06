source("Functions/mice_estimated_jags.R")

mice_estimated_jags_HB <- function(){
  mice.green <- mice_estimated_jags("Green Control")
  mice.henry <- mice_estimated_jags("Henry Control")
  mice.tea <- mice_estimated_jags("Tea Control")
  
  mice.mean <- matrix(NA, 3, length(mice.henry$mice.mean))
  mice.prec <- matrix(NA, 3, length(mice.henry$mice.prec))
  
  mice.mean[1, 1:length(mice.green$mice.mean)] <- mice.green$mice.mean
  mice.mean[2,] <- mice.henry$mice.mean
  mice.mean[3, 1:length(mice.tea$mice.mean)] <- mice.tea$mice.mean
  
  mice.prec[1, 1:length(mice.green$mice.prec)] <- mice.green$mice.prec
  mice.prec[2,] <- mice.henry$mice.prec
  mice.prec[3, 1:length(mice.tea$mice.prec)] <- mice.tea$mice.prec
  
  return <- list(mice.mean = mice.mean,
                 mice.prec = mice.prec)
  return(return)
}