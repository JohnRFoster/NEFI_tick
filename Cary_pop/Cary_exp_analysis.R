setwd("C:/Users/foste/Desktop/R_work/R/Cary")
library(tidyverse)
library(lubridate)

# Codes for "Status"
#   N <- new animal
#   P <- tagged in previous trapping session
#   RT <- retagged animal
#   R <- recaptured animal in current trapping session

# Codes for "Fate"
#   1 <- live animal of status P, RT, R
#   2 <- live animal of status N
#   3 <- dead animal of status P, RT, R
#   4 <- dead animal of status N (not tagged)
#   5 <- animal was relocated off trapping grid

mouse.df <- read.csv("Cary_mouse.csv", header = TRUE)
mouse.df$Full.Date.2 <- ymd(mouse.df$Full.Date.2)
mouse.df$Full.Date.1 <- ymd(mouse.df$Full.Date.1)
mouse.df$Tag.. <- as.character(mouse.df$Tag..)

all.mouse.df <- filter(mouse.df, Fate == 1 | Fate == 2)
new.mouse.df <- filter(mouse.df, Fate == 2)

count.all <- as.data.frame(table(all.mouse.df$Full.Date.1))
count.new <- as.data.frame(table(new.mouse.df$Full.Date.1))
day.of.year.all <- as.data.frame(table(all.mouse.df$day.of.year))


plot(count.all$Var1, count.all$Freq, main = "Total number of Mice Caught in Each Trapping Session")
plot(count.new$Var1, count.new$Freq, main = "New Mice Caught in Each Trapping Session")

plot(count.all$)

hist(count.all$Freq, main = "Frequency distribution of all mice caught per trapping session", 
     xlab = "Number caught per trapping session")
hist(count.new$Freq, main = "Frequency distribution of new mice caught per trapping session",
     xlab = "Number caught per trapping session")
