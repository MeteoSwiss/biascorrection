#' Bias dependent on initial condition
#' 
#' Computes mean de-biasing in dependence of initial conditions
#' 
##' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @details
#' This bias correction method assumes that the time-dependent mean bias depends on the
#' initial conditions. The method loosely follows the ideas outlined in Fuckar et al. (2014),
#' but in contrast to their approach, we use the forecast and observed conditions at the 
#' first day of the forecast as a proxy for the initial condition. Thereby, individual ensemble
#' members have varying bias correction depending on their respective initial conditions.
#' 
#' 
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
#' fcst.debias <- biascorrection:::initial(fcst[,1:20,], 
#' obs[,1:20], fcst.out=fcst, fc.time=fc.time, span=0.5)
#' 
#' 
#' @keywords util
initial <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.clim <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- sloess(obs.mn, span=span)
  oanom <- obs - obs.clim
  fcst.anom <- fcst - fcst.clim
  fcst.out.anom <- fcst.out - fcst.clim
  init <- oanom[1,] - mean(oanom[1,], na.rm=T)
  obs.fit <- oanom %*% init/ sum(init**2)
  ## all the ensemble members follow the same fit
  fcst.fit <- matrix(fcst.anom, nrow(fcst.anom)) %*% as.vector(fcst.anom[1,,]) / sum(fcst.anom[1,,]**2)
  ## predict the forecast anomalies with the model
  fcst.pred <- array((fcst.fit - obs.fit) %*% as.vector(fcst.out.anom[1,,]), dim(fcst.out.anom))
  ## piece the absolute debiased forecasts together
  fcst.debias <- fcst.out - (fcst.clim - obs.clim) - fcst.pred
  
  return(fcst.debias)
}
