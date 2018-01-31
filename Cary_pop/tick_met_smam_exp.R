library(tidyverse)
library(lubridate)
library(leaps)

t <- read.csv("tick_met")

###
lm.delta.doy <- lm(delta.nymph ~ day.of.year.y, data = t)
sum.lm.delta.doy <- summary(lm.delta.doy)
# Adj-R^2 = 0.06563

plot(x = t$day.of.year.y, y = t$delta.nymph)
abline(21.40299, -0.09896)

###
week.20 <- t %>% 
  filter(week == 20)

plot(x = week.20$gdd, y = week.20$delta.nymph)
