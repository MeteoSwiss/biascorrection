#' @name smooth_scale
#' 
#' @aliases smooth_scalespread
#' 
#' @title
#' Smoothed mean de-biasing with scaling
#' 
#' @description
#' Computes mean de-biasing with loess smoothing and adjusts variance. \code{smooth_scale}
#' scales the ensemble deviations from the climatology, whereas \code{smooth_scalespread}
#' only scales the ensemble anomalies (change to the ensemble spread).
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @details
#' The standard deviation (1/n definition) of the observations and simulations is computed against
#' the smoothed climatology  in order to get consistent results
#' should the forecast ensemble collapse on the observations. The disadvantage of this
#' approach is the dependence of the scaling on the climatological fit.
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(rnorm(215*30*51, mean=3, sd=0.2), c(215, 30, 51)) + 
#' 0.5*sin(seq(0,4,length=215))
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) + 
#' sin(seq(0,4, length=215))
#' fcst.debias <- biascorrection:::smooth_scale(fcst[,1:20,], obs[,1:20], fcst.out=fcst, span=0.5)
#' 
#' @keywords util
smooth_scale <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  fcst.clim <- sloess(fcst.mn, span=span)
  obs.clim <- sloess(obs.mn, span=span)
  obs.sd <- sqrt(apply((obs - obs.clim)**2, 1, mean, na.rm=T))
  obs.sdsmooth <- sloess(obs.sd, span=span)
  fcst.sd <- sqrt(apply((fcst - fcst.clim)**2, 1, mean, na.rm=T))
  fcst.sdsmooth <- sloess(fcst.sd, span=span)
  fcst.debias <- (fcst.out - fcst.clim) * obs.sdsmooth / fcst.sdsmooth + obs.clim
  return(fcst.debias)
}

#' @rdname smooth_scale
#' 
smooth_scalespread <-function (fcst, obs, fcst.out = fcst, 
                               span = min(1, 31/nrow(fcst)),  ...){
  fcst.ens <- rowMeans(fcst, dims = 2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims = 1, na.rm = T)
  obs.mn <- rowMeans(obs, dims = 1, na.rm = T)
  fcst.clim <- sloess(fcst.mn, span = span)
  obs.clim <- sloess(obs.mn, span = span)
  obs.sd <- sqrt(apply((obs - obs.clim)^2, 1, mean, na.rm = T))
  obs.sdsmooth <- sloess(obs.sd, span = span)
  fcst.sd <- sqrt(apply((fcst - fcst.clim)^2, 1, mean, na.rm = T))
  fcst.sdsmooth <- sloess(fcst.sd, span = span)
  fcst.out.ens <- rowMeans(fcst.out, dims=2, na.rm=T)
  fcst.debias <- (fcst.out - c(fcst.out.ens))* obs.sdsmooth / fcst.sdsmooth + 
    c(fcst.out.ens) - fcst.clim + obs.clim
  
  return(fcst.debias)
}
