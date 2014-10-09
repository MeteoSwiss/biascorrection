#' ccr
#' 
#' Computes climate conserving recalibration on smoothed monthly means
#' 
#' @param fcst array of forecast values (nyear, nlead, nens)
#' @param obs array of observations (nyear, nlead)
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @keywords util
#' @export
ccr <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  ## compute climatology
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  fcst.clim <- loess(fcst.mn ~ seq(along=fcst.mn), span=span)$fit
  obs.clim <- loess(obs.mn ~ seq(alon=obs.mn), span=span)$fit
  
  ## compute CCR prerequisites (Notation as in Weigel et al. 2009)
  x <- obs - obs.clim
  sig_x <- sqrt(apply(x**2, 1, mean))
  fi <- fcst - fcst.clim
  mu_f <- rowMeans(fi, dims=2)
  sig_mu <- sqrt(apply(mu_f**2, 1, mean))
  sig_ens <- sqrt(apply(apply(fi, 1:2, sd)**2, 1, mean))
  rho <- diag(cor(t(mu_f), t(x)))
  rr <- rho*sig_x / sig_mu
  ss <- sqrt(1 - rho**2) * sig_x / sig_ens
  ## put everything back together (with de-biasing)
  fi_out <- fcst.out - fcst.clim
  mu_fout <- rowMeans(fi_out, dims=2)
  fi_ccr <- as.vector(rr * mu_fout) + ss * (fi_out - as.vector(mu_fout)) + obs.clim  
  return(fi_ccr)
}
