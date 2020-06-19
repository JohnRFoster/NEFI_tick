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
  
  if(!is.null(met.proc)){
    if(met.proc == "vpd"){
      data$met <- data$vpd <- c(data$vpd[1,1], data$vpd[,2])
      data$met.mis <- which(is.na(data$met))
    } else if (met.proc == "min temp"){
      data$met <- c(data$temp.min[1,1], data$temp.min[,2])
      data$met.mis <- which(is.na(data$met))
    } else if (met.proc == "max temp"){
      data$met <- c(data$temp.max[1,1], data$temp.max[,2])
      data$met.mis <- which(is.na(data$met))
    } else if (met.proc == "min rh"){
      data$met <- c(data$rh.min[1,1], data$rh.min[,2])
      data$met.mis <- which(is.na(data$met))
    } else if (met.proc == "max rh"){
      data$met <- c(data$rh.max[1,1], data$rh.max[,2])
      data$met.mis <- which(is.na(data$met))
    } else if (met.proc == "preicp"){
      data$met <- c(data$precip[1,1], data$precip[,2])
      data$met.mis <- which(is.na(data$met))
    }  
  }
  
  data$temp.min <- c(data$temp.min[1,1], data$temp.min[,2])
  data$mis.temp.min <- which(is.na(data$temp.min))
  data$temp.min.range <- range(data$temp.min, na.rm = TRUE)
  data$met.mis.range <- range(data$met, na.rm = TRUE)
  
  data$temp.max <- data$rh.max <- data$rh.min <- data$precip <- data$vpd <- NULL
  data$mis.rh.max <- data$mis.rh.min <- data$mis.temp.max <- data$mis.vpd <- NULL
  
  return(data)
}

