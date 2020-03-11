S_E_quantiles <- function(mat, quants, percent = 1){
  store <- list()
  store$index.11 <- apply(mat[1,1,,], 1, quantile, quants)*percent  
  store$index.12 <- apply(mat[1,2,,], 1, quantile, quants)*percent  
  store$index.13 <- apply(mat[1,3,,], 1, quantile, quants)*percent
  store$index.21 <- apply(mat[2,1,,], 1, quantile, quants)*percent 
  store$index.22 <- apply(mat[2,2,,], 1, quantile, quants)*percent  
  store$index.23 <- apply(mat[2,3,,], 1, quantile, quants)*percent  
  store$index.31 <- apply(mat[3,1,,], 1, quantile, quants)*percent  
  store$index.32 <- apply(mat[3,2,,], 1, quantile, quants)*percent  
  store$index.33 <- apply(mat[3,3,,], 1, quantile, quants)*percent  
  return(store)
}
