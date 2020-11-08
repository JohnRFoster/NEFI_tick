library(plantecophys)
library(dplyr)
library(lubridate)

cary_weather_ensembles <- function(day.of.year){
  met <- read.csv("/projectnb/dietzelab/fosterj/Data/Cary_Met_Data_Daily.csv")
  met <- met %>% 
    mutate(DATE = mdy(DATE)) %>% 
    mutate(Year = year(DATE)) %>% 
    mutate(doy = yday(DATE)) %>% 
    mutate(vpd = RHtoVPD(MIN_RH, MIN_TEMP)) %>% 
    filter(DATE >= "1995-01-01") %>% 
    filter(DATE <= "2019-12-31") %>% 
    select(c("DATE","Year","doy","MIN_TEMP","MAX_TEMP","MAX_RH","MIN_RH","vpd"))
  
  calc_gdd <- function(base=10, max, min){
    gdd <- max(mean(max, min, na.rm = TRUE) - base, 0)
  }
  
  year <- met %>% 
    pull(Year) %>% 
    unique() 
  
  for(i in year){
    subset <- subset(met, Year == i)
    min.missing <- which(is.na(subset$MIN_TEMP))
    max.missing <- which(is.na(subset$MAX_TEMP))
    subset$MIN_TEMP[min.missing] <- mean(subset$MIN_TEMP, na.rm = TRUE)
    subset$MAX_TEMP[max.missing] <- mean(subset$MAX_TEMP, na.rm = TRUE)
    gdd <- rep(0, nrow(subset))
    for(t in 1:nrow(subset)){
      gdd[t] <- calc_gdd(10, subset$MAX_TEMP[t], subset$MIN_TEMP[t]) 
    }
    if(i == year[1]){
      cum.gdd <- cumsum(gdd)
    } else {
      cum.gdd <- c(cum.gdd, cumsum(gdd))
    }
  }
  
  met$cdd <- cum.gdd

  clim.ens <- list()
  for(ens in seq_along(year)){
    met.sample <- met %>%  
      select(-DATE) %>%
      filter(Year == year[ens]) %>% 
      filter(doy > day.of.year) %>% 
      filter(doy <= min(day.of.year + 35, 366)) %>% 
      rename(max.temp = MAX_TEMP,
             min.temp = MIN_TEMP,
             max.rh = MAX_RH,
             min.rh = MIN_RH,
             cum.gdd = cdd) %>% 
      select(-c(Year, doy))
    
    # some historical data is missing, set those days to the mean in the ensemble
    for(c in 1:ncol(met.sample)){
      if(any(is.na(met.sample[,c]))){
        fix <- which(is.na(met.sample[,c]))
        met.sample[fix,c] <- mean(met.sample[,c], na.rm = TRUE)
      }
    }
    clim.ens[[ens]] <- met.sample 
  }
  
  return(clim.ens)  
}


