get_month_effect <- function(t, m, ua, month.index){
  month.num <- month.index[t]
  
  if(any(grepl("alpha.month[1,", names(ua), fixed = TRUE))){ # if month.effect for each life stage
    larva.effect.name <- paste0("alpha.month[1,", month.num, "]")
    nymph.effect.name <- paste0("alpha.month[2,", month.num, "]")
    adult.effect.name <- paste0("alpha.month[3,", month.num, "]")
    month.effect <- c(ua[[larva.effect.name]][m], 
                      ua[[nymph.effect.name]][m],
                      ua[[adult.effect.name]][m])
    
  } else { # if month.effect for all life stages
    effect.name <- paste0("alpha.month[", month.num, "]")
    month.effect <- ua[[effect.name]][m]  
  }
  return(month.effect)
}