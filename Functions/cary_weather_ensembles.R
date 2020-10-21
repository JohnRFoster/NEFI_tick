library(plantecophys)
library(dplyr)
library(lubridate)

cary_weather_ensembles <- function(var){
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
    gdd <- max(mean(max, min) - base, 0)
  }
  
  year <- met %>% 
    pull(Year) %>% 
    unique() 
  
  n.years <- length(year)
  
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
  }
  
  met$cdd <- cum.gdd
  
  met.sample <- met %>%  
    select(all_of(c("Year", "doy", var))) %>% 
    pivot_wider(names_from = doy, 
                values_from = var)
  
}

