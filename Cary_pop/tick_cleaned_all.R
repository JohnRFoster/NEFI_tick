library(dplyr)
library(plyr)
library(lubridate)

## The below code is directly from tick_met_cleanup.R except here I am not adding the met data to 
## the tick data 

t <- read.csv("Cary_pop/CaryTick") 
t$Date <- as.character(t$Date)

t <- t[!grepl("1998", t$Date),] # delete 1998 data

t.1998 <- read.csv("Cary_pop/corrected 1998 tick density.csv")
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

tick <- tick[complete.cases(tick), ] # Rid tick data of NA values
tick$Date <- ymd(tick$Date)

tick <- union_all(tick, t.1998)
colnames(tick)[2] <- "DATE"

#########################
tick.green <- tick %>%
  filter(Grid == "Green Control")

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.green$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.green$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.green$n_adult), ncol = 1))

tick.green <- cbind(tick.green, delta.larvae, delta.nymph, delta.adult)
tick.green <- plyr::arrange(tick.green, DATE)

##
tick.green.e <- tick %>%
  filter(Grid == "Green Experimental")

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.green.e$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.green.e$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.green.e$n_adult), ncol = 1))

tick.green.e <- cbind(tick.green.e, delta.larvae, delta.nymph, delta.adult)
tick.green.e <- plyr::arrange(tick.green.e, DATE)
##
tick.henry <- tick %>%
  filter(Grid == "Henry Control")


delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.henry$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.henry$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.henry$n_adult), ncol = 1))

tick.henry <- cbind(tick.henry, delta.larvae, delta.nymph, delta.adult)
tick.henry <- plyr::arrange(tick.henry, DATE)
##
tick.henry.e <- tick %>%
  filter(Grid == "Henry Experimental")

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.henry.e$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.henry.e$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.henry.e$n_adult), ncol = 1))

tick.henry.e <- cbind(tick.henry.e, delta.larvae, delta.nymph, delta.adult)
tick.henry.e <- plyr::arrange(tick.henry.e, DATE)
##
tick.tea <- tick %>%
  filter(Grid == "Tea Control") 

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.tea$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.tea$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.tea$n_adult), ncol = 1))

tick.tea <- cbind(tick.tea, delta.larvae, delta.nymph, delta.adult)
tick.tea <- plyr::arrange(tick.tea, DATE)
##
tick.tea.e <- tick %>%
  filter(Grid == "Tea Experimental")

delta <- as.matrix(0, ncol = 1)
delta.larvae <- rbind(delta, as.matrix(diff(tick.tea.e$n_larvae), ncol = 1))
delta.nymph <- rbind(delta, as.matrix(diff(tick.tea.e$n_nymph), ncol = 1))
delta.adult <- rbind(delta, as.matrix(diff(tick.tea.e$n_adult), ncol = 1))

tick.tea.e <- cbind(tick.tea.e, delta.larvae, delta.nymph, delta.adult)
tick.tea.e <- plyr::arrange(tick.tea.e, DATE)
##########################


tick.all <- rbind(tick.green, tick.green.e, tick.henry, tick.henry.e, tick.tea, tick.tea.e)
write.csv(tick.all, file = "Cary_pop/tick_cleaned")

