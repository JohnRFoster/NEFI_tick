## Function to grab survival posteriors to use as priors in tick population models

#' @param larva.driver met driver for larva survival
#' @param nymph.driver met driver for nymph survival

get_survival <- function(larva.driver, nymph.driver){
  
  # get larva
  if(is.null(larva.driver)){
    load("../FinalOut/SurvivalModels/Flat_Larva/Flat_Larvae_NULL_AllTicks_NoRandom.RData")
    larva.survival.mean <- mean(as.matrix(jags.out))
    larva.survival.prec <- 1 / var(as.matrix(jags.out))
    larva.beta.mean <- NULL
    larva.beta.prec <- NULL
    larva.beta.diff.mean <- NULL
    larva.beta.diff.prec <- NULL
  } else if (larva.driver == "max temp"){
    load("../FinalOut/SurvivalModels/Flat_Larva/Flat_Larvae_temp_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    larva.survival.mean <- mean(out[,"beta[1]"])
    larva.survival.prec <- 1 / var(out[,"beta[1]"])
    larva.beta.mean <- mean(out[,"beta[2]"])
    larva.beta.prec <- 1 / var(out[,"beta[2]"])
    larva.beta.diff.mean <- mean(out[,"beta[3]"])
    larva.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else if (larva.driver == "max rh"){
    load("../FinalOut/SurvivalModels/Flat_Larva/Flat_Larvae_rh_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    larva.survival.mean <- mean(out[,"beta[1]"])
    larva.survival.prec <- 1 / var(out[,"beta[1]"])
    larva.beta.mean <- mean(out[,"beta[2]"])
    larva.beta.prec <- 1 / var(out[,"beta[2]"])
    larva.beta.diff.mean <- mean(out[,"beta[3]"])
    larva.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else if (larva.driver == "vpd") {
    load("../FinalOut/SurvivalModels/Flat_Larva/Flat_Larvae_vpd_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    larva.survival.mean <- mean(out[,"beta[1]"])
    larva.survival.prec <- 1 / var(out[,"beta[1]"])
    larva.beta.mean <- mean(out[,"beta[2]"])
    larva.beta.prec <- 1 / var(out[,"beta[2]"])
    larva.beta.diff.mean <- mean(out[,"beta[3]"])
    larva.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else {
    stop("Did not set larva.driver! Must be NULL, temp, rh, or vpd")
  }
  
  # get nymph
  if(is.null(nymph.driver)){
    load("../FinalOut/SurvivalModels/Flat_Nymph/Flat_Nymph_NULL_AllTicks_NoRandom.RData")
    nymph.survival.mean <- mean(as.matrix(jags.out))
    nymph.survival.prec <- 1 / var(as.matrix(jags.out))
    nymph.beta.mean <- NULL
    nymph.beta.prec <- NULL
    nymph.beta.diff.mean <- NULL
    nymph.beta.diff.prec <- NULL
  } else if (nymph.driver == "max temp") {
    load("../FinalOut/SurvivalModels/Flat_Nymph/Flat_Nymph_temp_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    nymph.survival.mean <- mean(out[,"beta[1]"])
    nymph.survival.prec <- 1 / var(out[,"beta[1]"])
    nymph.beta.mean <- mean(out[,"beta[2]"])
    nymph.beta.prec <- 1 / var(out[,"beta[2]"])
    nymph.beta.diff.mean <- mean(out[,"beta[3]"])
    nymph.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else if (nymph.driver == "max rh") {
    load("../FinalOut/SurvivalModels/Flat_Nymph/Flat_Nymph_rh_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    nymph.survival.mean <- mean(out[,"beta[1]"])
    nymph.survival.prec <- 1 / var(out[,"beta[1]"])
    nymph.beta.mean <- mean(out[,"beta[2]"])
    nymph.beta.prec <- 1 / var(out[,"beta[2]"])
    nymph.beta.diff.mean <- mean(out[,"beta[3]"])
    nymph.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else if (nymph.driver == "vpd") {
    load("../FinalOut/SurvivalModels/Flat_Nymph/Flat_Nymph_vpd_Above_AllTicks_NoRandom.RData")
    out <- as.matrix(jags.out)
    nymph.survival.mean <- mean(out[,"beta[1]"])
    nymph.survival.prec <- 1 / var(out[,"beta[1]"])
    nymph.beta.mean <- mean(out[,"beta[2]"])
    nymph.beta.prec <- 1 / var(out[,"beta[2]"])
    nymph.beta.diff.mean <- mean(out[,"beta[3]"])
    nymph.beta.diff.prec <- 1 / var(out[,"beta[3]"])
  } else {
    stop("Did not set nymph.driver! Must be NULL, temp, rh, or vpd")
  }

  survival.post <- list(larva.survival.mean = larva.survival.mean,
                        larva.survival.prec = larva.survival.prec,
                        larva.beta.mean = larva.beta.mean,
                        larva.beta.prec = larva.beta.prec,
                        larva.beta.diff.mean = larva.beta.diff.mean,
                        larva.beta.diff.prec = larva.beta.diff.prec,
                        nymph.survival.mean = nymph.survival.mean,
                        nymph.survival.prec = nymph.survival.prec,
                        nymph.beta.mean = nymph.beta.mean,
                        nymph.beta.prec = nymph.beta.prec,
                        nymph.beta.diff.mean = nymph.beta.diff.mean,
                        nymph.beta.diff.prec = nymph.beta.diff.prec)
  return(survival.post)
}
