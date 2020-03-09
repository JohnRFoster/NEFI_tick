source("Functions/known_states.R")
source("Functions/ch_cary.R")

mouse_data_jags <- function(site){
  tea <- suppressWarnings(ch_cary(site))
  aug <- matrix(0, nrow(tea), ncol(tea))
  y <- rbind(tea, aug)
  
  tea.ks <- known_states(tea)
  x <- rbind(tea.ks, aug)
  
  cap.date <- colnames(tea)
  cap.date <- as.Date(cap.date,format="%Y-%m-%d")
  start.date <- cap.date[1]
  end.date <- cap.date[length(cap.date)]
  
  diff <- diff.Date(cap.date)
  dt <- as.numeric(c(diff,0))
  dt.index <- cumsum(dt)
  
  met <- read.csv("../Met_Cary")
  met.start <- which(met$DATE == as.character(start.date))
  met.end <- which(met$DATE == as.character(end.date))
  met <- met[met.start:met.end,]
  
  precip <- met[,"TOT_PREC"]
  temp <- met[,"MAX_TEMP"]
  temp <- scale(temp, scale = FALSE)
  temp.mis <- which(is.na(temp))
  
  data <- list(y = y,  
               dt = dt.index, 
               ind = nrow(y), 
               time = ncol(y),
               precip = as.matrix(precip),
               temp = as.matrix(temp),
               temp.mis = temp.mis,
               days = sum(dt))
  
  return(data)
}
