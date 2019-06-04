library(ggplot2)
library(grid)
library(gridExtra)

#' This function plots each chain from JAGS as it's own distribution
#'
#' @param posterior mcmc list object with iterations from JAGS
#' @param params character vector of parameters (or states) to plot. Default to all.
#' @export
#' @examples plot_posterior_chain(out$params)

plot_posterior_chains <- function(posterior,params=colnames(posterior[[1]])){
  mcmc <- as.matrix(posterior)
  n.chains <- length(posterior)
  chain <- factor(rep(1:n.chains,each=nrow(mcmc)/n.chains))
  p <- list()
  mcmc <- as.matrix(mcmc[,params])
    for(i in 1:ncol(mcmc)){
      data <- data.frame(chain = chain, param = mcmc[,i])
      title <- colnames(mcmc)[i]
      p[[i]] <- ggplot(data, aes(x = param,fill = chain)) + 
        geom_density(alpha = 0.3) +
        labs(title = title)
  }
  do.call(grid.arrange,p)
}