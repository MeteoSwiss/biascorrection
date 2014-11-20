#' qqmap_semi
#' 
#' Computes bias correction with semi-parametric quantile mapping
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of qqmaping (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @details
#' The quantile-quantile mapping implemented here is based on the assumption that
#' the observation data are approximately normally distributed. In contrast, quantiles
#' for the forecasts are directly estimated based on the assumption that there are enough
#' degrees of freedom to estimate the quantiles directly. This may strongly affect forecasts
#' in regions with skill. 
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(runif(215*30*51), c(215, 30, 51)) + 
#' 0.5*sin(seq(0,4,length=215)) 
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) + 
#' sin(seq(0,4, length=215)) 
#' fcst.debias <- biascorrection:::qqmap_semi(fcst[,1:20,], obs[,1:20], fcst.out=fcst, span=0.5)
#' 
#' @keywords util
qqmap_semi <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  ## compute climatology of obs
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- sloess(obs.mn, span=span)
  ## compute standard deviation
  obs.sd <- sqrt(rowMeans((obs - obs.clim)**2, dims=1, na.rm=T))
  obs.sdsmooth <- sloess(obs.sd, span=span)
  ## compute quantile-quantile mapping
  ## fcst.quant <- pnorm(fcst.out, mean=fcst.clim, sd=fcst.sdsmooth)
  fcst.rank <- array(t(apply(fcst.out, 1, rank)), dim(fcst.out))
  fcst.quant <- (fcst.rank - 0.5) / apply(fcst.rank, 1, max)
  fcst.debias <- qnorm(fcst.quant, mean=obs.clim, sd=obs.sdsmooth)  
  return(fcst.debias)
}
