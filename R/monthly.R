#' Daily calibration with monthly mean
#' 
#' Computes monthly mean de-biasing
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for monthly aggregation
#' @param fcout.time forecast time for array to which bias correction is applied
#' for back compatibility with leave-one-out cross-validation in \code{\link{debias}}
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @keywords util
monthly <- function(fcst, obs, fcst.out=fcst, fc.time, fcout.time=fc.time, ...){
  if (length(fcout.time) != length(fcst.out[,,1])) {
    stop('Time (fcout.time) is not of correct dimension/length')
  }
  if (length(fc.time) < length(fcst[,,1])){
    stop('Not enough forecast time steps (fc.time) supplied')
  }
  fcst.ens <- rowMeans(fcst, dims=2)
  ## compute monthly bias
  bias.mon <- rowMeans(month(fcst.ens - obs, fc.time), dims=1)
  ## de-bias with monthly means
  fcst.debias <- fcst.out - as.vector(bias.mon[format(fcout.time, '%m')])
  return(fcst.debias)
}
