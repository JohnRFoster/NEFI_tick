library(tidyverse)
library(lubridate)


mice_06_18 <- function(grid){
  mice.data.06_08 <- read.csv("../Mouse_mast_raw_trapping_data_2006_2018.csv", 
                          header = TRUE,
                          na.strings=c(c("", " ", "      "),"NA"),
                          skip = 2)
  
  smam.df <- mice.data.06_08
  smam.df$Date <- as.character(smam.df$Date)
  smam.df$Year <- as.character(smam.df$Year)
  
  # smam.df <- tidyr::separate(smam.df, Date, sep = "[[:alpha:]]", c("month1", "day"))
  
  smam.df$Date <- gsub("Sept", "Sep", smam.df$Date) 
  smam.df$Date <- gsub("SEPT", "Sep", smam.df$Date)
  smam.df$Date <- gsub("September", "Sep", smam.df$Date)
  smam.df$Date <- gsub("Sepember", "Sep", smam.df$Date)
  
  
  smam.df <- smam.df %>% unite(Year, Date, sep = " ", col = "Full.Date.2", remove = FALSE)
  
  
  date.formats <- c("%Y %d-%d %b", "%Y %d-%d %B", "%Y %d - %d %B", "%Y %d - %d %b", "%Y %d/%m/%Y", "%Y %d-%d %b.",
                    "%Y %d-%d%B", "%Y %d&%d %B", "%Y %d&%d %b", "%Y %d-%b", "%d-%B", "%Y %d,%d %b", "%Y %d %b- $d %b.",
                    "%Y %d %b-%d %b", "%Y %d-%d %B/ %b", "%Y %d,%d %B", "%Y %b %d-%d", "%Y %B %d-%d", "%Y %b. %d-%d",
                    "%Y %b.%d-%d", "%Y %d %b-%b %d", "%Y %B %d-%b %d", "%Y %b %d-%b %d", "%Y %b %d-%b %d",
                    "%Y %B %d-%b %d", "%Y %b. %d-%d", "%Y %b. %d", "%Y %b %d -%d", "%Y %B%d-%d", 
                    "%Y %m/%d-%d/%Y", "%Y %m/%d/-%m/%d/%Y", "%Y %B %d-","%Y %B %d","%Y %m/%d, %d/%Y","%Y %B/%d-%d/%Y")
  
  all.dates.func <- function(data, formats){
    a <- list()
    for(i in 1:length(formats)){
      a[[i]] <- as.Date(data, format = formats[i])
      a[[1]][!is.na(a[[i]])] <- a[[i]][!is.na(a[[i]])]
    }
    a[[1]]
  }
  
  smam.df$Full.Date.2 <- all.dates.func(smam.df$Full.Date.2, date.formats)
  # length(which(is.na(smam.df$Full.Date.2)))
  na.index <- which(is.na(smam.df$Full.Date.2))
  smam.df$Full.Date.2 <- as.character(smam.df$Full.Date.2)
  smam.df$Full.Date.2[na.index] <- "2010/10/01"
  
  bad.index <- which(smam.df$Date == "8/31-9/1/2006")
  smam.df$Full.Date.2[bad.index] <- "2006/09/01"
  
  bad.index <- which(smam.df$Date == "10/31-11/1/2006")
  smam.df$Full.Date.2[bad.index] <- "2006/11/01"
  
  bad.index <- which(smam.df$Date == "10/16-17/07")
  smam.df$Full.Date.2[bad.index] <- "2007/10/17"
  
  bad.index <- which(smam.df$Date == "11/6-7/07")
  smam.df$Full.Date.2[bad.index] <- "2007/11/07"
  
  bad.index <- which(smam.df$Date == "7/31-8/1/2007")
  smam.df$Full.Date.2[bad.index] <- "2007/08/01"
  
  bad.index <- which(smam.df$Date == "10/30-31/07")
  smam.df$Full.Date.2[bad.index] <- "2007/10/31"
  
  bad.index <- which(smam.df$Date == "10/18-19/07")
  smam.df$Full.Date.2[bad.index] <- "2007/10/19"
  
  bad.index <- which(smam.df$Date == "7/31-8/1/2008")
  smam.df$Full.Date.2[bad.index] <- "2008/08/01"
  
  bad.index <- which(smam.df$Date == "9/30-10/1/2008")
  smam.df$Full.Date.2[bad.index] <- "2008/10/01"
  
  bad.index <- which(smam.df$Date == "7/27-28/10")
  smam.df$Full.Date.2[bad.index] <- "2010/07/28"
  
  bad.index <- which(smam.df$Date == "8/31-9/1/2010")
  smam.df$Full.Date.2[bad.index] <- "2010/09/01"
  
  bad.index <- which(smam.df$Date == "6/23-6/24/2011")
  smam.df$Full.Date.2[bad.index] <- "2011/06/24"
  
  bad.index <- which(smam.df$Date == "5/31-6/1/2011")
  smam.df$Full.Date.2[bad.index] <- "2011/06/01"
  
  bad.index <- which(smam.df$Date == "6/30-7/1/2011")
  smam.df$Full.Date.2[bad.index] <- "2011/07/01"
  
  bad.index <- which(smam.df$Date == "10/31-11/1/2012")
  smam.df$Full.Date.2[bad.index] <- "2012/11/01"
  
  bad.index <- which(smam.df$Date == "5/31-6/1/2012")
  smam.df$Full.Date.2[bad.index] <- "2012/06/01"
  
  bad.index <- which(smam.df$Date == "7/31-8/1/2012")
  smam.df$Full.Date.2[bad.index] <- "2012/08/01"
  
  bad.index <- which(smam.df$Date == "10/31-11/1/2013")
  smam.df$Full.Date.2[bad.index] <- "2013/11/01"
  
  bad.index <- which(smam.df$Date == "6/13-6/14/2013")
  smam.df$Full.Date.2[bad.index] <- "2013/06/14"
  
  bad.index <- which(smam.df$Date == "7/1-7/2/2013")
  smam.df$Full.Date.2[bad.index] <- "2013/07/02"
  
  
  # bad <- filter(smam.df, Full.Date.2 < as.Date("2005-01-01"))
  # head(bad)

  smam.df$Full.Date.2 <- ymd(smam.df$Full.Date.2)
  smam.df <- smam.df %>% mutate(Full.Date.1 = Full.Date.2 - 1) %>% mutate(day.of.year = yday(Full.Date.1))
  
  smam <- smam.df[,c("Grid",
                     "Full.Date.1",
                     "Full.Date.2",
                     "Day.1",
                     "Day.2",
                     "Tag..",
                     "Fate")]
  
  alive <- c(1, 2)                                    # codes for alive individuals
  #smam <- smam %>% filter(Fate %in% alive)            # extract only live individuals
  smam <- subset(smam, Fate == 1 | Fate == 2)
  
  smam[4:5] <- as.integer(!is.na(smam[4:5]))          # converts trap histories to 1 or 0
  
  smam$Full.Date.1 <- as.Date(smam$Full.Date.1)           # convert factors to dates
  smam$Full.Date.2 <- as.Date(smam$Full.Date.2)           # convert factors to dates
  
  m <- subset(smam, Grid == grid)
  m <- m[,c("Tag..", "Day.1", "Day.2", "Full.Date.1", "Full.Date.2")]
  m <- subset(m, !is.na(Tag..))
  
  m$Tag.. <- as.character(m$Tag..)
  m.ls <- split(m, m$Full.Date.1)                     # split on first capture date for sampling occasion
  for(i in 1:length(m.ls)){
    m.ls[[i]] <- m.ls[[i]][c(-4, -5)]
  }
  
  day.1 <- as.character(unique(m$Full.Date.1))        # 1st capture date of sampling occasion
  day.2 <- as.character(unique(m$Full.Date.2))        # 2nd capture date of sampling occasion
  days <- c(rbind(day.1, day.2))                      # vector of unique trapping days (for colnames)
  
  ch.base <- merge(m.ls[[1]], m.ls[[2]], by = "Tag..", all = TRUE)
  l.m <- length(m.ls)*1
  for(i in 3:l.m){                                    # loop through the rest
    g <- as.data.frame(m.ls[[i]])
    ch.base <- merge(ch.base, g, by = "Tag..", all = TRUE)
  }
  ch.base <- as.matrix(ch.base[,-1])                            # convert all NAs to 0
  for(i in 1:nrow(ch.base)){
    for(t in 1:ncol(ch.base)){
      if(is.na(ch.base[i,t])){
        ch.base[i,t] <-0 
        }
      ch.base[i,t] <- as.numeric(ch.base[i,t])
    }
  }
  colnames(ch.base) <- days
  
  # known_states <- function(ch){
  #   state <- ch
  #   for (i in 1:dim(ch)[1]){
  #     n1 <- min(which(ch[i,] != 0))
  #     n2 <- max(which(ch[i,] != 0))
  #     if(n2 > n1){state[i, n1:n2] <- 1}
  #     #state[i, n1:n2] <- 1
  #   }
  #   return(state)
  # }
  # 
  # ks <- known_states(ch.base)
  # return <- apply(ks, 2, sum)
  
  return(list(full.matrix = ch.base,
              table = m))
}
