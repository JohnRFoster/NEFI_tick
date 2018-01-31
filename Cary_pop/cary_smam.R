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


data <- read.csv("Mouse_mast_raw_trapping_data.csv", header = TRUE, sep = ",") # Mouse trapping data all sites all years
smam.df <- data
smam.df$Date <- as.character(smam.df$Date)
smam.df$Year <- as.character(smam.df$Year)


smam.df$Date <- gsub("Sept", "Sep", smam.df$Date) 
smam.df$Date <- gsub("SEPT", "Sep", smam.df$Date)
smam.df$Date <- gsub("September", "Sep", smam.df$Date)
smam.df$Date <- gsub("Sepember", "Sep", smam.df$Date)
smam.df$Date <- gsub("July 26-37", "July 26-27", smam.df$Date)

smam.df <- smam.df %>% unite(Year, Date, sep = " ", col = "Full.Date.2", remove = FALSE)

date.formats <- c("%Y %d-%d %b", "%Y %d-%d %B", "%Y %d - %d %B", "%Y %d - %d %b", "%Y %d/%m/%Y", "%Y %d-%d %b.",
                  "%Y %d-%d%B", "%Y %d&%d %B", "%Y %d&%d %b", "%Y %d-%b", "%d-%B", "%Y %d,%d %b", "%Y %d %b- $d %b.",
                  "%Y %d %b-%d %b", "%Y %d-%d %B/ %b", "%Y %d,%d %B", "%Y %b %d-%d", "%Y %B %d-%d", "%Y %b. %d-%d",
                  "%Y %b.%d-%d", "%Y %d %b-%b %d", "%Y %B %d-%b %d", "%Y %b %d-%b %d", "%Y %b %d-%b %d",
                  "%Y %B %d-%b %d", "%Y %b. %d-%d", "%Y %b. %d", "%Y %b %d -%d", "%Y %B%d-%d")

all.dates.func <- function(data, formats){
  a <- list()
  for(i in 1:length(formats)){
    a[[i]] <- as.Date(data, format = formats[i])
    a[[1]][!is.na(a[[i]])] <- a[[i]][!is.na(a[[i]])]
  }
  a[[1]]
}

smam.df$Full.Date.2 <- all.dates.func(smam.df$Full.Date.2, date.formats)
smam.df$Full.Date.2[23080] <- as.Date("2004-10-28")
smam.df$Full.Date.2[23081] <- as.Date("2004-11-01")
#which(is.na(smam.df$Full.Date.2))

smam.df$Full.Date.2 <- ymd(smam.df$Full.Date.2)
smam.df <- smam.df %>% mutate(Full.Date.1 = Full.Date.2 - 1) %>% mutate(day.of.year = yday(Full.Date.1))





write.csv(smam.df, "Cary_mouse.csv", sep = ",")

green.ctrl <- filter(smam.df, Grid == "Green Control")           # the following six filter() calls are making 
green.exp <- filter(smam.df, Grid == "Green Experimental")       # data frames for each site, all years
henry.ctrl <- filter(smam.df, Grid == "Henry Control")
henry.exp <- filter(smam.df, Grid == "Henry Experimental")
tea.ctrl <- filter(smam.df, Grid == "Tea Control")
tea.exp <- filter(smam.df, Grid == "Tea Experimental")
