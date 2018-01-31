library(tidyverse)
library(lubridate)
library(leaps)

m <- read.csv("Met_Cary")

GrowDegDay <- function(temp, Tbase){
  
  # temp = vector of temperatures to be used in growing degree day calculation
  # Tbase = temperature threshold for growing degree days
  # cum = function, True = cumulative, False = single calculation
  
  for (i in 1:length(temp)) {
    temp[is.na(temp)] <- 0
    if(temp[i] > Tbase){
      gdd[i] <- temp[i] - Tbase
    } else {
      gdd[i] <- 0
    }
  }
  gdd <- cumsum(gdd)
  return(gdd)
}
   
     


met <- m
met$DATE <- ymd(met$DATE)

met <- met %>% 
          mutate(day.of.year = yday(DATE)) %>% 
          mutate(year = year(DATE)) %>%
          mutate(week = week(DATE))

y.1 <- met %>% 
  filter(year == 1995) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.2 <- met %>% 
  filter(year == 1996) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.3 <- met %>% 
  filter(year == 1997) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.4 <- met %>% 
  filter(year == 1998) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.5 <- met %>% 
  filter(year == 1999) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.6 <- met %>% 
  filter(year == 2000) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.7 <- met %>% 
  filter(year == 2001) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.8 <- met %>% 
  filter(year == 2002) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.9 <- met %>% 
  filter(year == 2003) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.10 <- met %>% 
  filter(year == 2004) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))
y.11 <- met %>% 
  filter(year == 2005) %>% 
  select(AVE_TEMP) %>%
  mutate(gdd = GrowDegDay(AVE_TEMP, 10))

cum.gdd <- bind_rows(y.1, y.2, y.3, y.4, y.5, y.6, y.7, y.8, y.9, y.10, y.11)
met <- cbind(met, cum.gdd)


t <- read.csv("CaryTick")

t.1998 <- read.csv("corrected 1998 tick density.csv")
t.1998$Date <- t.1998$Date %>% 
  as.character() %>% 
  as.Date("%m/%d/%Y") %>% 
  ymd()

t.1998 <- t.1998 %>% mutate(day.of.year = yday(t.1998$Date))

t.1998$X <- as.character(t.1998$X)

for(i in 1:length(t.1998$X)){
  if(t.1998$X[i] == "Green Con"){t.1998$X[i] <- "Green Control"}
  if(t.1998$X[i] == "Green Exp"){t.1998$X[i] <- "Green Experimental"}
  if(t.1998$X[i] == "Henry Con"){t.1998$X[i] <- "Henry Control"}
  if(t.1998$X[i] == "Henry Exp"){t.1998$X[i] <- "Henry Experimental"}
  if(t.1998$X[i] == "Tea Con"){t.1998$X[i] <- "Tea Control"}
  if(t.1998$X[i] == "Tea Exp"){t.1998$X[i] <- "Tea Experimental"}
}

colnames(t.1998)[1] <- "Grid"
colnames(t.1998)[3:5] <- c("Larvae.m2", "Nymphs.m2", "Adults.m2")

tick <- t
# compute absolute number of ticks caught by multiplying density by area (450 m^2)of tick drags
tick <- tick %>%                          
          select(-X, -X.1, -Year) %>%
          mutate(n_larvae = round(Larvae.m2 * 450, 0)) %>%
          mutate(n_nymphs = round(Nymphs.m2 * 450, 0)) %>%
          mutate(n_adults = round(Adults.m2 * 450, 0))

tick <- tick[complete.cases(tick), ]
tick$Date <- ymd(tick$Date)

tick <- union_all(tick, t.1998)
colnames(tick)[2] <- "DATE"

#########################
tick.green <- tick %>%
  filter(Grid == "Green Control") %>%
  inner_join(., met, by = "DATE") %>%
  select(-day.of.year.x, -X) 

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.green$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.green$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.green$n_adult), ncol = 1))

tick.green <- cbind(tick.green, delta.larvae, delta.nymph, delta.adult)
tick.green <- tick.green[complete.cases(tick.green), ]

tick.henry <- tick %>%
  filter(Grid == "Henry Control") %>% 
  inner_join(., met, by = "DATE") %>%
  select(-day.of.year.x, -X) 

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.henry$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.henry$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.henry$n_adult), ncol = 1))

tick.henry <- cbind(tick.henry, delta.larvae, delta.nymph, delta.adult)
tick.henry <- tick.henry[complete.cases(tick.henry), ]

tick.tea <- tick %>%
  filter(Grid == "Tea Control") %>%
  inner_join(., met, by = "DATE") %>%
  select(-day.of.year.x, -X) 
  
delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.tea$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.tea$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.tea$n_adult), ncol = 1))

tick.tea <- cbind(tick.tea, delta.larvae, delta.nymph, delta.adult)
tick.tea <- tick.tea[complete.cases(tick.tea), ]
##########################


tick.all <- rbind(tick.green, tick.henry, tick.tea)
write.csv(tick.all, file = "tick_met")


