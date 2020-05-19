library(lubridate)

future_met <- function(site, gdd.base = 10){
  met <- read.csv("../Cary_Met_Data_Daily.csv")
  met <- met[, c("DATE","MAX_TEMP", "MIN_TEMP", "MAX_RH","TOT_PREC")]
  met$DATE <- as.Date(as.character(met$DATE), format = c("%m/%d/%Y"))
  
  sites <- c("Green Control","Henry Control","Tea Control")
  N_site <- length(sites)               # number of sites
  raw.dat <- read.csv("../tick_cleaned")   # read in data
  raw.dat$DATE <- as.character(raw.dat$DATE) # convert to date
  
  # find last day in timeseries for each site
  last.day <- vector()
  for(i in 1:length(sites)){
    subset <- subset(raw.dat, Grid == sites[i])
    last.day[i] <- as.character(subset$DATE[nrow(subset)])
  }
  
  # assign 
  if (site == "Green Control") s <- 1
  if (site == "Henry Control") s <- 2
  if (site == "Tea Control") s <- 3
  
  # index for last day for site being forecasted
  last.index <- which(met$DATE == last.day[s])
  met.seq <- last.index:which(met$DATE == "2017-12-31") # all 2018 data is NA
  
  met <- met[met.seq,]
  
  met$DATE <- as.Date(met$DATE)
  year <- format(as.Date(met$DATE, format="%Y-%m-%d"),"%Y")[met.seq]
  year <- as.numeric(as.factor(year))
  
  met$year <- format(as.Date(met$DATE, format="%Y-%m-%d"),"%Y")
  met$year <- as.numeric(as.factor(met$year))

  calc_gdd <- function(base, max, min){
    gdd <- max(mean(max, min) - base, 0)
  }
  
  gdd.vec <- vector()
  cum.gdd <- vector()
  for(i in 1:length(unique(met$year))){
    subset <- subset(met, year == i)
    min.missing <- which(is.na(subset$MIN_TEMP))
    max.missing <- which(is.na(subset$MAX_TEMP))
    subset$MIN_TEMP[min.missing] <- mean(subset$MIN_TEMP, na.rm = TRUE)
    subset$MAX_TEMP[max.missing] <- mean(subset$MAX_TEMP, na.rm = TRUE)
    gdd <- vector()
    for(t in 1:nrow(subset)){
      gdd[t] <- calc_gdd(gdd.base, subset$MAX_TEMP[t], subset$MIN_TEMP[t]) 
    }
    cum.gdd <- c(cum.gdd, cumsum(gdd))
    gdd.vec <- c(gdd.vec, gdd)
  }
  
  # if na's set to mean
  temp.scale <- scale(met$MAX_TEMP, scale = FALSE)
  temp <- met$MAX_TEMP
  if(any(is.na(temp))){
    na.temp <- which(is.na(temp))
    temp.scale[na.temp] <- mean(temp.scale, na.rm = TRUE)
    temp[na.temp] <- mean(temp, na.rm = TRUE)
  }
  precip <- met$TOT_PREC
  if(any(is.na(precip))){
    na.precip <- which(is.na(precip))
    precip[na.precip] <- mean(precip, na.rm = TRUE)
  }
  rh.scale <- scale(met$MAX_RH, scale = FALSE)
  rh <- met$MAX_RH
  if(any(is.na(rh))){
    na.rh <- which(is.na(rh))
    rh.scale[na.rh] <- mean(rh.scale, na.rm = TRUE)
    rh[na.rh] <- mean(rh, na.rm = TRUE)
  }
  
  # create data frame
  future.met <- data.frame(temp = temp,
                           temp.scale = temp.scale,
                           precip = precip,
                           rh = rh,
                           rh.scale = rh.scale,
                           cum.gdd = cum.gdd,
                           year.index = met$year,
                           date = met$DATE)
  
  return(future.met)
}

# test.green <- future_met("Green Control")
# test.henry <- future_met("Henry Control")
# test.tea <- future_met("Tea Control")
