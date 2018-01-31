setwd("C:/Users/foste/Desktop/R_work/R/Cary")
library(tidyverse)
library(sp)
library(geosphere)
library(raster)
library(lubridate)

# Green control polygon
coords.green <- matrix(c(-73.767730, 41.793760,
                         -73.769390, 41.793530, 
                         -73.769980, 41.794560, 
                         -73.768243, 41.795100,
                         -73.767730, 41.793760), 
                     ncol = 2, byrow = TRUE)

# Henry control polygon
coords.henry <- matrix(c(-73.721642, 41.803156,
                         -73.721813, 41.804577,
                         -73.720461, 41.804711,
                         -73.720048, 41.803354, 
                         -73.721642, 41.803156), 
                       ncol = 2, byrow = TRUE)

# Henry control polygon
coords.tea <- matrix(c(-73.732226, 41.793532,
                       -73.731523, 41.794803,
                       -73.730064, 41.794025,
                       -73.730772, 41.793065, 
                       -73.732226, 41.793532), 
                       ncol = 2, byrow = TRUE)

area.calc <- function(coords){
  p <- Polygon(coords)
  sp <- SpatialPolygons(list(Polygons(list(p), ID = "a")),
                        proj4string = CRS("+init=epsg:4269"))
  area <- areaPolygon(sp)
  return(area)
}

# area is m^2
area.green <- area.calc(coords.green)  
area.henry <- area.calc(coords.henry) 
area.tea <- area.calc(coords.tea) 


## normaloze by mean and standard deviation 
## look at three sites independently first then work into heirarchy 






# mice data, .all = all live captures, .new = only new mice caught
data.mice <- read.csv("Cary_mouse.csv")

green.m.all <- filter(data.mice, Grid == "Green Control") %>% filter(Fate == 1 | Fate == 2)
green.m.new <- filter(data.mice, Grid == "Green Control") %>% filter(Fate == 2)
green.all <- as.data.frame(table(green.m.all$Full.Date.1)) %>% 
                filter(Freq != 0) %>%
                mutate(density = Freq / area.green * 1000) %>%
                mutate(std_density = density / sd(density))

green.new <- as.data.frame(table(green.m.new$Full.Date.1)) %>%
                filter(Freq != 0) %>%
                mutate(density = Freq / area.green * 1000) %>%
                mutate(std_density = density / sd(density))

henry.m.all <- filter(data.mice, Grid == "Henry Control") %>% filter(Fate == 1 | Fate == 2)
henry.m.new <- filter(data.mice, Grid == "Henry Control") %>% filter(Fate == 2)
henry.all <- as.data.frame(table(henry.m.all$Full.Date.1)) %>% 
                filter(Freq != 0) %>%
                mutate(density = Freq / area.henry)
henry.new <- as.data.frame(table(henry.m.new$Full.Date.1)) %>% 
                filter(Freq != 0) %>%
                mutate(density = Freq / area.henry)

tea.m.all <- filter(data.mice, Grid == "Tea Control") %>% filter(Fate == 1 | Fate == 2)
tea.m.new <- filter(data.mice, Grid == "Tea Control") %>% filter(Fate == 2)
tea.all <- as.data.frame(table(tea.m.all$Full.Date.1)) %>% 
                filter(Freq != 0) %>%
                mutate(density = Freq / area.tea)
tea.new <- as.data.frame(table(tea.m.new$Full.Date.1)) %>% 
                filter(Freq != 0) %>%
                mutate(density = Freq / area.tea)

# tick data intake
data <- read.csv("Mouse_mast_Tick_Drag.csv", header = TRUE)
tick.df <- data

tick.df$Date <- as.character(tick.df$Date)
tick.df[135, 3] <- "6/12/02"                  # these three rows have two sampling dates (6/12/02-6/13/02)
tick.df[290, 3] <- "6/23/03"                  # coerced to just first day of the sampling effort 
tick.df[445, 3] <- "10/25/04"
tick.df$Date <- mdy(tick.df$Date)
tick.df <- tick.df %>% mutate(day.of.year = yday(Date)) # add day-of-year column to data frame

# data frames for each control site
green.nymph <- tick.df %>% 
  filter(Grid == "Green Control") %>% 
  dplyr::select(Date, Nymphs.m2) %>% 
  mutate(std_density = Nymphs.m2 / sd(Nymphs.m2))
  

henry.nymph <- tick.df %>% 
  filter(Grid == "Henry Control") %>% 
  select(Date, Nymphs.m2) %>% 
  mutate(per_1000m2 = Nymphs.m2 * 100)

tea.nymph <- tick.df %>% 
  filter(Grid == "Tea Control") %>% 
  select(Date, Nymphs.m2) %>% 
  mutate(per_1000m2 = Nymphs.m2 * 100)   

## data visualization


ggplot() +
  geom_line(data = green.all, aes(x = as.Date(Var1), y = std_density), color = "blue") +
  geom_line(data = green.nymph, aes(x = as.Date(Date), y = std_density), color = "red")

ggplot() +
  geom_line(data = green.new, aes(x = as.Date(Var1), y = std_density), color = "blue") +
  geom_line(data = green.nymph, aes(x = as.Date(Date), y = std_density), color = "red")









### the following was to check for a solid polygon, and to see all sites spatially
# p.green <- Polygon(coords.green)
# sp.green <- SpatialPolygons(list(Polygons(list(p.green), ID = "a")),
                            # proj4string = CRS("+proj=longlat +datum=WGS84"))
# plot(sp.green, axes = TRUE)

# p.henry <- Polygon(coords.henry)
# sp.henry <- SpatialPolygons(list(Polygons(list(p.henry), ID = "a")),
                            # proj4string = CRS("+proj=longlat +datum=WGS84"))
# plot(sp.henry, axes = TRUE)

# p.tea <- Polygon(coords.tea)
# sp.tea <- SpatialPolygons(list(Polygons(list(p.tea), ID = "a")),
                            # proj4string = CRS("+proj=longlat +datum=WGS84"))
# plot(sp.tea, axes = TRUE)

# spl <- list(sp.green, sp.tea, sp.henry)
# m <- do.call(bind, spl)
# plot(m, axes = TRUE)
###