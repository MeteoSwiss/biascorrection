#' unbias
#' 
#' Computes simple mean de-biasing
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @keywords util
unbias <- function(fcst, obs, fcst.out=fcst, ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  fcst.debias <- fcst.out - fcst.mn + obs.mn
  
  return(fcst.debias)
}
