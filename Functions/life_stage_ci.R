life.stage.ci <- function(pred, type, quants = c(0.025, 0.5, 0.975)){
  if(type == "forecast"){
    time <- 1:(length(pred[1,,1]))
  } else {
    time <- 1:(length(pred[1,,1])-1)
  }
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