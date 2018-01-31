library(tidyverse)
library(mra)
library(lubridate)

data <- read.csv("cary_mouse.csv")

# subset mouse data to just what is needed for the mark-recapture analysis
# specifically: individual, date, capture histories, and site (grid)
mr <- data 
mr$Status <- as.integer(mr$Status)
dead <- c(3, 6, 11, 16)
mr <- mr %>% 
  filter(!Status %in% dead) %>%     # delete captures of dead mice
  select(Tag.., Full.Date.1, Day.1, Day.2, Grid) 

# convert capture histories to 1 (captured) or 0 (not captured)
mr$Day.1 <- as.numeric(mr$Day.1)
for (i in 1:length(mr$Day.1)){
  if(mr$Day.1[i] == 1){
    mr$Day.1[i] <- 0
  } else {
    mr$Day.1[i] <- 1
  }
}
mr$Day.2 <- as.numeric(mr$Day.2)
for (i in 1:length(mr$Day.2)){
  if(mr$Day.2[i] == 1){
    mr$Day.2[i] <- 0
  } else {
    mr$Day.2[i] <- 1
  }
}

# deleting captures that do not have a Tag (mouse is unidentified)
mr$Tag.. <- as.character(mr$Tag..)
mr <- mr %>% 
  filter(Tag.. != "")





