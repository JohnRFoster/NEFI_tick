library(tidyverse)
library(lubridate)

dat <- read.csv("Cary_pop/tick_cleaned") # tick data
gr <- read.csv("GreenTicks.csv", header = FALSE) # green control mice data

green <- dat %>% select(c(Grid, DATE, n_larvae, n_nymphs, n_adults)) %>% 
  filter(Grid == "Green Control") # get only tick data from green control
green$DATE <- ymd(green$DATE)

mice <- read.csv("Cary_pop/Cary_mouse.csv") %>% 
  filter(Grid == "Green Control")
mice$Full.Date.1 <- ymd(mice$Full.Date.1)

# get year-month for each sampling day (ticks and mice)
day.mice <- unique(mice$Full.Date.1) # unique sampling days: mice
day.tick <- unique(green$DATE) # unique sampling days: tick

year.t1 <- vector()
for(t in 1:length(day.tick)){
  up <- day.tick[t] - 365 + 15
  down <- day.tick[t] - 365 - 15
  for(m in 1:length(day.mice)){
    if((day.mice[m] >= down) && (day.mice[m] <= up)){
      year.t1[t] <- as.character(day.mice[m])      
    } 
  }
}

year.t2 <- vector()
for(t in 1:length(day.tick)){
  up <- day.tick[t] - 365*2 + 15
  down <- day.tick[t] - 365*2 - 15
  for(m in 1:length(day.mice)){
    if((day.mice[m] >= down) && (day.mice[m] <= up)){
      year.t2[t] <- as.character(day.mice[m])      
    } 
  }
}
cbind(as.character(day.tick), year.t1, year.t2)

