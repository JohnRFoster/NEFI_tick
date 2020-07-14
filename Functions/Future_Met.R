library(lubridate)

future_met <- function(site, gdd.base = 10){
  met <- read.csv("../Cary_Met_Data_Daily.csv")
  met <- met[, c("DATE","MAX_TEMP", "MIN_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")]
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
  
  met <- met[1:which(met$DATE == "2017-12-31"),]
  
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
  
  # if na's set to mean and scale, return both
  scale_mean <- function(vec){
    scale.vec <- scale(vec, scale = FALSE)
    if(any(is.na(vec))){
      na <- which(is.na(vec))
      scale.vec[na] <- mean(scale.vec, na.rm = TRUE)
      vec[na] <- mean(vec, na.rm = TRUE)
    }
    return(list(var.scale = as.vector(scale.vec),
                var = vec))
  }
  max.temp <- scale_mean(met$MAX_TEMP)
  min.temp <- scale_mean(met$MIN_TEMP)
  max.rh <- scale_mean(met$MAX_RH)
  min.rh <- scale_mean(met$MIN_RH)
  precip <- scale_mean(met$TOT_PREC)

  # create data frame
  future.met <- data.frame(max.temp = max.temp$var[met.seq],
                           max.temp.scale = max.temp$var.scale[met.seq],
                           min.temp = min.temp$var[met.seq],
                           min.temp.scale = min.temp$var.scale[met.seq],
                           max.rh = max.rh$var[met.seq],
                           max.rh.scale = max.rh$var.scale[met.seq],
                           min.rh = min.rh$var[met.seq],
                           min.rh.scale = min.rh$var.scale[met.seq],
                           precip = precip$var[met.seq],
                           cum.gdd = cum.gdd[met.seq],
                           year.index = met$year[met.seq],
                           date = met$DATE[met.seq])
  
  return(future.met)
}

# test.green <- future_met("Green Control")
# test.henry <- future_met("Henry Control")
# test.tea <- future_met("Tea Control")
