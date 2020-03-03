library(plantecophys)
site_data_met <- function(site, met.variable, data, obs.model = TRUE){
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
    if(met.variable == "max temp"){
      dat.jags$met <- data$temp.max[,s]
      dat.jags$met.mis <- data$mis.temp.max[,s]
    } else if(met.variable == "min temp"){
      dat.jags$met <- data$temp.min[,s]
      dat.jags$met.mis <- data$mis.temp.min[,s]
    } else if(met.variable == "vpd"){
      dat.jags$met <- data$vpd[,s]
      dat.jags$met.mis <- data$mis.vpd[,s]
    } else if(met.variable == "max rh"){
      dat.jags$met <- data$rh.max[,s]
      dat.jags$met.mis <- data$mis.rh.max[,s]
    } else if(met.variable == "min rh"){
      dat.jags$met <- data$rh.min[,s]
      dat.jags$met.mis <- data$mis.rh.min[,s]
    } else if(met.variable == "precip"){
      dat.jags$met <- data$precip[,s]
    }
    dat.jags$met.range <- range(dat.jags$met, na.rm = TRUE)
  }
  
  
  
  if(obs.model){
    obs.index <- c(1, dat.jags$dt.index)
    met.obs <- data$temp.min[obs.index[1:dat.jags$N_est],s]
    met.obs.miss <- which(is.na(met.obs))
    met.obs.range <- range(met.obs, na.rm = TRUE)
    
    dat.jags$met.obs <- met.obs
    dat.jags$met.obs.miss <- met.obs.miss
    dat.jags$met.obs.range <- met.obs.range
  }
  return(dat.jags)
}