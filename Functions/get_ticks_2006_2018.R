#' This function retrieves tick data from 2006-2018
#'
#' @param site the site for which ticks are needed
#' @export
#' @examples get_ticks_2006_2018("Green Control")


library(tidyverse)
library(xlsx)

get_ticks_2006_2018 <- function(site){
  t <- read.xlsx("../Mous_mast_Tick_Drag_2006_2018_dateEdit.xlsx", sheetIndex = 1, startRow = 3) 
  
  # compute absolute number of ticks caught by multiplying density by area (450 m^2)of tick drags
  individual.ticks <- t[,-7] %>% # 7th column is just NA
    dplyr::filter(Grid %in% site) %>%  # filter to site
    dplyr::mutate(n_larvae = round(Larvae.m2 * 450, 0)) %>% # calculate ticks caught
    dplyr::mutate(n_nymphs = round(Nymphs.m2 * 450, 0)) %>%
    dplyr::mutate(n_adults = round(Adults.m2 * 450, 0)) %>% 
    dplyr::select(c("Date", "n_larvae", "n_nymphs", "n_adults")) 
  
  df <- diff.Date(individual.ticks$Date)
  
  return <- list(individual.ticks = t(individual.ticks[,-1]),
                 df = df)
  
  return(return)
}
