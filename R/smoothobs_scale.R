#' smoothobs_scale
#' 
#' Computes mean de-biasing with loess smoothing (obs. only) and adjusts variance
#' 
#' @param fcst array of forecast values (nyear, nlead, nens)
#' @param obs array of observations (nyear, nlead)
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(rnorm(215*30*51, mean=3, sd=0.2), c(215, 30, 51)) + 
#' 0.5*sin(seq(0,4,length=215))
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) + 
#' sin(seq(0,4, length=215))
#' fcst.debias <- smoothobs_scale(fcst[,1:20,], obs[,1:20], fcst.out=fcst, span=0.5)
#' 
#' @keywords util
#' @export
smoothobs_scale <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- sloess(obs.mn, span=span)
  obs.sd <- apply(obs, 1, sd, na.rm=T)
  obs.sdsmooth <- sloess(obs.sd, span=span)
  fcst.sd <- apply(fcst, 1, sd, na.rm=T)
  fcst.debias <- (fcst.out - fcst.mn) * obs.sdsmooth / fcst.sd + obs.clim
  return(fcst.debias)
}
