new_mouse_observation <- function(capt.history, mice, days, index) {
  current.date <- days[index]
  
  if (current.date %in% mice$Full.Date.1) {
    current.matrix <- mice %>%
      filter(Full.Date.1 == current.date) %>%
      select(Tag.., Day.1) %>%
      filter(Day.1 == 1)
  } else if (current.date %in% mice$Full.Date.2) {
    current.matrix <- mice %>%
      filter(Full.Date.2 == current.date) %>%
      select(Tag.., Day.2)
  }
  
  # full capture history
  capt.history <-
    full_join(capt.history, current.matrix, by = "Tag..")
  
  # convert NAs to 0
  last <- ncol(capt.history)
  capt.history[which(is.na(capt.history[, last])), last] <- 0
  
  # capt.history <- capt.history[,c(1,ncol(capt.history))]
  
  return(capt.history = capt.history)
}