hb_restructure <- function(data, met.proc){
  
  # restructure cdd for only one vector
  data$gdd <- c(data$gdd[1,1], data$gdd[2,])
  
  # site one (Green Control) is one day ahead of the other two
  # sites, so we need to add one to keep dt.index consistent with
  # gdd and met
  data$dt.index[2:3,] <- data$dt.index[2:3,] + 1
  
  # sequence for permuting matrices for each site
  seq.days <- array(NA, 
                    dim = c(max(data$N_est-1), # max number of days to estimate
                            max(data$df, na.rm = TRUE), # max days between sampling events
                            length(data$N_est))) # number of sites
  
  for(s in 1:length(data$N_est)){
    for(i in 1:(data$N_est[s]-1)){
      xx <- (data$dt.index[s,i+1]-1):data$dt.index[s,i]
      seq.days[i,1:length(xx),s] <- xx
    }
  }    
  data$seq.days <- seq.days  
  
  data$temp.min <- c(data$temp.min[1,1], data$temp.min[,2])
  data$temp.max <- c(data$temp.max[1,1], data$temp.max[,2])
  data$rh.min <- c(data$rh.min[1,1], data$rh.min[,2])
  data$rh.max <- c(data$rh.max[1,1], data$rh.max[,2])
  data$precip <- c(data$precip[1,1], data$precip[,2])
  data$vpd <- c(data$vpd[1,1], data$vpd[,2])
  
  data$mis.rh.max <- which(is.na(data$rh.max))
  data$mis.rh.min <- which(is.na(data$rh.min))
  data$mis.temp.max <- which(is.na(data$temp.max))
  data$mis.temp.min <- which(is.na(data$temp.min))
  data$mis.vpd <- which(is.na(data$vpd))
  
  return(data)
}

