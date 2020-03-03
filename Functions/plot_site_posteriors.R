library(ggplot2)
library(grid)
library(gridExtra)

#' This function plots each chain from JAGS as it's own distribution
#'
#' @param all.params matrix (or data frame) of site posteriors with all chains
#' @param n.site the number of sites within all.params, default = 3
#' @export
#' @examples plot_posterior_chain(out$params)
#' 

plot_site_posteriors <- function(all.params, n.site = 3){
  
  all.params <- as.data.frame(all.params)
  
  n.iter <- nrow(all.params)/n.site
  site <- factor(rep(1:n.site, each = n.iter))
  
  p <- list()
  for(i in 1:ncol(all.params)){
    data <- data.frame(site = site, param = all.params[,i])
    title <- colnames(all.params)[i]
    p[[i]] <- ggplot(data, aes(x = param, fill = site)) + 
      geom_density(alpha = 0.3) +
      labs(title = title) 
  }
  do.call(grid.arrange,p)
}