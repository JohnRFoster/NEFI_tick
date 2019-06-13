
site_data_met <- function(site, met.variable, data){
  # index for site
  dat.jags <- list()
  if(site == "Green Control"){
    s <- 1
    dat.jags$y <- data$y[,-73,1]
  } else if(site == "Henry Control"){
    s <- 2
    dat.jags$y <- data$y[,,2]
  } else {
    s <- 3
    dat.jags$y <- data$y[,-73,3]
  }
  
  dat.jags$N_est <- data$N_est[s]
  dat.jags$N_days <- data$N_days[s]
  dat.jags$dt.index <- data$dt.index[s,]
  dat.jags$df <- data$df[s,]
  dat.jags$gdd <- data$gdd[s,]
  # dat.jags$year <- data$year[,s]
  
  if(!is.null(met.variable)){
    if(met.variable == "temp"){
      dat.jags$met <- data$met[,1,s]
      dat.jags$temp.mis <- data$temp.mis[,s]
    } else if(met.variable == "rh"){
      dat.jags$met <- data$met[,2,s]
      dat.jags$rh.mis <- data$rh.mis[,s]
    } else {
      dat.jags$met <- data$met[,3,s]
    }
  }
  return(dat.jags)
}