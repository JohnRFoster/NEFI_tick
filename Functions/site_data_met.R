library(plantecophys)
library(dplyr)
library(lubridate)

site_data_met <- function(site, met.variable, data, dir = "", 
                          obs.model = TRUE, time.effect = NULL){
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
  
  if(!is.null(time.effect)){
    raw.dat <- read.csv("/projectnb/dietzelab/fosterj/Data/tick_cleaned")
    raw.dat$DATE <- as.Date(raw.dat$DATE) # convert to date
    
    dates <- raw.dat %>% 
      filter(Grid == site) %>% 
      select(DATE)
    if(time.effect == "month"){
      dat.jags$month.index <- lubridate::month(dates[,1])
      dat.jags$n.months <- unique(lubridate::month(dates[,1]))
    }  
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
  
  seq.days <- matrix(NA, dat.jags$N_est - 1, max(dat.jags$df, na.rm = TRUE))
  for (i in 1:(dat.jags$N_est - 1)) {
    xx <- (dat.jags$dt.index[i + 1] - 1):dat.jags$dt.index[i]
    seq.days[i, 1:length(xx)] <- xx
  }
  dat.jags$seq.days <- seq.days
  
  dat.jags$R <- diag(1, 3, 3)
  
  return(dat.jags)
}