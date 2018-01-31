library(tidyverse)
library(mra)
library(lubridate)

data <- read.csv("Cary_pop/cary_mouse.csv")

# subset mouse data to just what is needed for the mark-recapture analysis
# specifically: individual, date, capture histories, and site (grid)
mr <- data 
mr$Status <- as.integer(mr$Status)
dead <- c(3, 6, 11, 16)
mr <- mr %>% 
  select(Tag.., Full.Date.1, Day.1, Day.2, Grid, Status, Sex) 

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

mr.1 <- mr %>% 
  mutate(Bin = 1:nrow(mr))

# change all day 1 dead captures to 2
day.1 <- mr.1 %>% 
  filter(Day.1 == 1 & Day.2 == 0 & Status %in% dead) %>% 
  mutate(Day.1, Day.1 = 2) 

mr.1[match(day.1$Bin, mr.1$Bin), ] <- day.1 # replace fixed capture history into mr.1

# change all day 2 dead captures to 2
day.2 <- mr.1 %>% 
  filter(Day.1 == 0 & Day.2 == 1 & Status %in% dead) %>% 
  mutate(Day.2, Day.2 = 2)

mr.1[match(day.2$Bin, mr.1$Bin), ] <- day.2 # replace fixed capture history into mr.1

# change all mice captured on day 1 and dead on day 2 to 2
both <- mr.1 %>% 
  filter(Day.1 == 1 & Day.2 == 1 & Status %in% dead) %>% 
  mutate(Day.2, Day.2 = 2)

mr.1[match(both$Bin, mr.1$Bin), ] <- both # replace fixed capture history into mr.1


# these 3 rows have capture histories of 0, 0; but also dead. They all have ear tags, so I made their 
# day 1 capture as a 2
mr.1[7897, 3] <- 2
mr.1[23080, 3] <- 2
mr.1[23081, 3] <- 2

# deleting captures that do not have a Tag (mouse is unidentified)
mr.1$Tag.. <- as.character(mr.1$Tag..)
mr <- mr %>% 
  filter(Tag.. != "")

history.1 <- mr.1 %>% 
  filter(Full.Date.1 == "1995-04-25") %>% 
  select(Day.1, Day.2) %>% 
  as.matrix()

xy <- F.cjs.covars(nrow(history.1), ncol(history.1))
for(j in 1:ncol(history.1)){ assign(paste("x", j, sep = ""), xy$x[,,j]) }
dipper.cjs <- F.cjs.estim( ~x1+x2, ~x1+x2, history.1 )
dipper.cjs
plot(dipper.cjs)





h2 <- mr.1 %>% 
  filter(Full.Date.1 == "1995-05-23") %>% 
  select(Day.1, Day.2) %>% 
  as.matrix()

xy <- F.cjs.covars(nrow(h2), ncol(h2))
for(j in 1:ncol(h2)){ assign(paste("x", j, sep = ""), xy$x[,,j]) }
h2.cjs <- F.cjs.estim( ~x1+x2, ~x1+x2, h2 )
h2.cjs
plot.cjs(h2.cjs)
