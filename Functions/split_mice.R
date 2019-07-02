#' This splits out$params into the parameters and mice estimates
#'
#'@param out model output 
#'

split_mice <- function(out){
  model <- list(params = NULL, mice = NULL)
  mfit <- as.matrix(out$params, chains = TRUE)
  mice <- grep("mice[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  model$mice <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, mice)])
  model$params <- ecoforecastR::mat2mcmc.list(mfit[, -mice])
  return(model)
}
