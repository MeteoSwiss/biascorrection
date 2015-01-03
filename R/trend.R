#' Bias with linear time trend
#' 
#' Computes mean de-biasing with linear time trend
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for trend aggregation
#' @param fcout.time forecast time for array to which bias correction is applied
#' for back compatibility with leave-one-out cross-validation in \code{\link{debias}}
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @details
#' This bias correction method assumes that the bias can be decomposed 
#' into a stationary seasonal cycle (as in method \code{\link{smoothobs}}) 
#' and a linear time trend estimated from the residuals.
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(rnorm(215*30*51), c(215, 30, 51)) + 
#' 0.5*sin(seq(0,4,length=215)) + 
#' rep(seq(0,1,length=30), each=215)
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) + 
#' sin(seq(0,4, length=215)) + 
#' rep(seq(0,3,length=30), each=215)
#' fc.time <- outer(1:215, 1981:2010, function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x)
#' fcst.debias <- biascorrection:::trend(fcst[,1:20,], 
#' obs[,1:20], fcst.out=fcst, fc.time=fc.time, span=0.5)
#' 
#' 
#' @keywords util
trend <- function(fcst, obs, fcst.out=fcst, fc.time, fcout.time=fc.time, span=min(1, 31/nrow(fcst)), ...){
  if (length(fcout.time) != length(fcst.out[,,1])) {
    stop('Time (fcout.time) is not of correct dimension/length')
  }
  if (length(fc.time) < length(fcst[,,1])){
    stop('Not enough forecast time steps (fc.time) supplied')
  }
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  fcst.clim <- sloess(fcst.mn, span=span)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- sloess(obs.mn, span=span)
  ## bias <- fcst.clim - obs.clim
  
  ## compute smoothing for bias in both lead time and fcst year
  ## anom <- fcst.ens - obs - bias
  anom <- fcst.ens - obs - fcst.clim + obs.clim
  lead <- rep(1:nrow(anom), ncol(anom))
  year <- rep(as.numeric(format(fc.time[1,seq(1,ncol(fcst))], '%Y')), each=nrow(fcst))
  anom.lm <- lm(as.vector(anom) ~ poly(lead, 2)*year)
  
  ## compute debiased forecast
  leadout <- rep(1:nrow(fcst.out), ncol(fcst.out))
  yearout <- rep(as.numeric(format(fcout.time[1,], '%Y')), each=nrow(fcst.out))
  fcst.debias <- fcst.out - 
    predict(anom.lm, newdata=list(lead=leadout, year=yearout)) +
    obs.clim - fcst.clim

  return(fcst.debias)
}
