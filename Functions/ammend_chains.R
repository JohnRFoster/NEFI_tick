#' This takes the continued jags run and ammends it to the previous run
#'
#'@param dir file path to model output
#'@param num.chains integer, the number of chains run in parallel
#'@param num.out integer, the number of succesive files written out (the minimum between chains)
#'@param thin number of iterations to keep, default 5000
#'@param save do you want to save the combined chains? defult is no (NULL) 
#'

source('Functions/combine_chains.R')

ammend_chains <- function(dir, num.chains, num.out, thin = 5000, save = FALSE){
  
  all.chains <- list()
  for(c in 1:num.chains){
    c1 <- list()
    for(i in 1:num.out){
      load <- paste(dir, "_c", c, "_", i, ".RData", sep = "")
      load(load)
      if(num.out == 1){
        c1$params <- as.mcmc(out$out$params)
        c1$predict <- as.mcmc(out$out$predict)
        c1$m.cols <- as.mcmc(out$out$m.cols)
      } else {
        c1$params <- rbind(c1$params, as.mcmc(out$out$params))
        c1$predict <- rbind(c1$predict, as.mcmc(out$out$predict))
        c1$m.cols <- rbind(c1$m.cols, as.mcmc(out$out$m.cols))
      }
    }
    iter <- dim(c1$params)[1]
    if(iter <= 5000){
      all.chains$params[[c]] <- as.mcmc(c1$params)
      all.chains$predict[[c]] <- as.mcmc(c1$predict)
      all.chains$m.cols[[c]] <- as.mcmc(c1$m.cols)
    } else {
      iter <- round(seq(1, by = iter/thin, length = thin))
      all.chains$params[[c]] <- as.mcmc(c1$params[iter,])
      all.chains$predict[[c]] <- as.mcmc(c1$predict[iter,])
      all.chains$m.cols[[c]] <- as.mcmc(c1$m.cols[iter,])
    }
  }
  
  out.test <- list(params = as.mcmc(all.chains$params),
                   predict = as.mcmc(all.chains$predict),
                   m.cols = as.mcmc(all.chains$m.cols))
  
  if(is.character(save)){
    save(out.test, file = paste(dir, "_AllChains.RData", sep = ""))
  }
  return(out.test)
}


# load("../FinalOut/HB_Partial_GDD/Temp/WindowLoop/GDDSwitch_Window_Temp_c1_1.RData")
# load("../FinalOut/HB_Partial_GDD/Temp/WindowLoop/GDDSwitch_Window_Temp_c1_2.RData")
# c1 <- list()
# c1[[1]] <- out
# c1[[2]] <- out
# 
# num.chains <- 5
# num.out <- 3
# dir <- "../FinalOut/HB_Partial_GDD/Temp/WindowLoop"
# model <- "GDDSwitch_Window_Temp"
# 
# part1 <- combine_chains(file.path(dir,model))









