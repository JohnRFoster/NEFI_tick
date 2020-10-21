library(plantecophys)
library(dplyr)
scale_met_forecast <- function(){
  met <- read.csv("/projectnb/dietzelab/fosterj/Data/Cary_Met_Data_Daily.csv")
  met$DATE <- lubridate::mdy(met$DATE)
  met <- met %>% 
    mutate(vpd = RHtoVPD(met$MIN_RH, met$MIN_TEMP)) %>% 
    filter(DATE >= "1995-01-01") %>% 
    filter(DATE <= "2005-12-31") %>% 
    select(c("MIN_TEMP","MAX_TEMP","MAX_RH","MIN_RH","vpd")) 
  
  met.means <- apply(met, 2, mean, na.rm = TRUE)
  
  return(met.means)
} 