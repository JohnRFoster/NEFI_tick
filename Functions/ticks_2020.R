library(lubridate)
library(tidyverse)

ticks_2020 <- function(site){
  
  dat <- read.csv("/projectnb/dietzelab/fosterj/Data/Tick_drag_May-Sept_2020.csv", 
                  header = TRUE, skip = 1)
  
  # replace Date column name
  dat.col <- grep("Date", colnames(dat))
  colnames(dat)[dat.col] <- "Date"
  
  dat$Date <- ymd(dat$Date) # convert to date
  dat$Grid <- as.character(dat$Grid) # convert to character
  
  # remove "X" column and subset grid
  dat <- dat %>% 
    select(c(Date, Grid, Larvae, Nymphs, Adult.total)) %>% 
    filter_all(all_vars(!is.na(.))) %>% 
    filter(Grid == site) 
  
  return(dat)
} 