#' This function plots posterior and prior distributions
#'
#' @param posterior vector of posterior mcmc iterations
#' @param dist named distribution used as prior - like you would call in R
#' @param args named list of arguments for dist
#' @export
#' @examples plot_post_prior(params[,"k.l2n.low"], dnorm, list(mean = 400, sd = 1/sqrt(0.0001)))

library(ggplot2)

plot_post_prior <- function(posterior, dist, args){
  data <- data.frame(posterior = posterior)
  ggplot <- ggplot(data, aes(x = posterior)) +
    stat_function(fun = dist, 
                  n = nrow(data), 
                  args = args,
                  geom = "area",
                  color = "black",
                  alpha = 0.3,
                  aes(fill = "Prior")) +
    geom_density(aes(x = posterior, fill = "Posterior"), alpha = 0.3) +
    scale_fill_manual(name = "",
                      values = c("blue", "red")) +
    labs(x = "Value", y = "") +
    theme_classic()
  return(ggplot)
}
