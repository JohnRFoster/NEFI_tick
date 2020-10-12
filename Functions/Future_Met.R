library(lubridate)

future_met <- function(site, gdd.base = 10){
  met <- read.csv("/projectnb/dietzelab/fosterj/Data/Cary_Met_Data_Daily.csv")
  met <- met[, c("DATE","MAX_TEMP", "MIN_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")]
  met$DATE <- as.Date(as.character(met$DATE), format = c("%m/%d/%Y"))
  
  ## cumulative gdd calculation
  met$year <- year(met$DATE)
  
  year <- met %>% 
    pull(year) %>% 
    unique()
  
  calc_gdd <- function(base=10, max, min){
    gdd <- max(mean(max, min) - base, 0)
  }
  
  gdd.vec <- vector()
  cum.gdd <- vector()
  for(i in year){
    subset <- subset(met, year == i)
    min.missing <- which(is.na(subset$MIN_TEMP))
    max.missing <- which(is.na(subset$MAX_TEMP))
    subset$MIN_TEMP[min.missing] <- mean(subset$MIN_TEMP, na.rm = TRUE)
    subset$MAX_TEMP[max.missing] <- mean(subset$MAX_TEMP, na.rm = TRUE)
    gdd <- vector()
    for(t in 1:nrow(subset)){
      gdd[t] <- calc_gdd(10, subset$MAX_TEMP[t], subset$MIN_TEMP[t]) 
    }
    cum.gdd <- c(cum.gdd, cumsum(gdd))
    gdd.vec <- c(gdd.vec, gdd)
  }
  
  met$cdd <- gdd.vec
  
  raw.dat <- read.csv("/projectnb/dietzelab/fosterj/Data/tick_cleaned")   # read in tick data
  raw.dat$DATE <- as.character(raw.dat$DATE) # convert to date
  
  raw.dat <- filter(raw.dat, Grid == site)
  last.day <- raw.dat$DATE[nrow(raw.dat)]
  
  met <- met %>% 
    filter(DATE >= last.day)

  return(met)
}

# test.green <- future_met("Green Control")
# test.henry <- future_met("Henry Control")
# test.tea <- future_met("Tea Control")
