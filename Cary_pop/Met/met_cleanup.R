library(tidyverse)
library(lubridate)

temp <- dir(pattern = "*.csv")
data <- do.call(rbind,lapply(temp, read.csv))

years <- c("1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005")

var <- c("DATE", "MIN_TEMP", "MAX_TEMP", "AVE_TEMP", "MIN_RH", "MAX_RH", "AVE_RH", "TOT_PREC")

m <- select(data, var)
met <- select(data, var)
met$DATE <- mdy(met$DATE)

met <- met %>% 
  mutate(day.of.year = yday(DATE)) %>% 
  mutate(year = year(DATE)) %>%
  mutate(week = week(DATE))

# no max or min 1539:1541, 2342 (no data at all), 3099. All other rows that have both max and min, so 
# average is calculated as the mean between the max and min temps. 

for (i in 1:length(met$AVE_TEMP)) {
  if(is.na(met$AVE_TEMP[i])){
    met$AVE_TEMP[i] <- apply(met[i, 2:3], 1, mean)
  }
}

for (i in 1:length(met$AVE_RH)) {
  if(is.na(met$AVE_RH[i])){
    met$AVE_RH[i] <- apply(met[i, 5:6], 1, mean)
  }
}

write.csv(met, file = "Met_Cary")
