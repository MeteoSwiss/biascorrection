#' ccr_monthly
#' 
#' Computes climate conserving recalibration on smoothed monthly means
#' 
#' @param fcst array of forecast values (nyear, nlead, nens)
#' @param obs array of observations (nyear, nlead)
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for monthly aggregation
#' @param fcout.time forecast time for array to which bias correction is applied
#' for back compatibility with leave-one-out cross-validation in \code{\link{debias}}
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @keywords util
#' @export
ccr_monthly <- function(fcst, obs, fcst.out=fcst, fc.time, fcout.time=fc.time, ...){
  if (length(fcout.time) != length(fcst.out[,,1])) {
    stop('Time (fcout.time) is not of correct dimension/length')
  }
  if (length(fc.time) < length(fcst[,,1])){
    stop('Not enough forecast time steps (fc.time) supplied')
  }
  ## compute monthly forecasts and observations
  fcst.mon <- month(fcst, fc.time[seq(along=fcst[,,1])])
  obs.mon <- month(obs, fc.time[seq(obs)])
  
  ## compute climatologies
  fcst.clim <- rowMeans(fcst.mon, dims=1, na.rm=T)
  obs.clim <- rowMeans(obs.mon, dims=1, na.rm=T)
  
  ## compute ccr_monthly prerequisites (Notation as in Weigel et al. 2009)
  x <- obs.mon - obs.clim
  sig_x <- sqrt(apply(x**2, 1, mean))
  fi <- fcst.mon - fcst.clim
  mu_f <- rowMeans(fi, dims=2)
  sig_mu <- sqrt(apply(mu_f**2, 1, mean))
  sig_ens <- sqrt(apply(apply(fi, 1:2, sd)**2, 1, mean))
  rho <- diag(cor(t(mu_f), t(x)))
  rr <- rho*sig_x / sig_mu
  ss <- sqrt(1 - rho**2) * sig_x / sig_ens
  
  ## put everything back together
  monstr.out <- format(fcout.time, '%m')
  fi_out <- fcst.out - fcst.clim[monstr.out]
  mu_fout <- rowMeans(fi_out, dims=2)
  fi_ccr <- as.vector(rr[monstr.out]*mu_fout) + 
    ss[monstr.out]*(fi_out - as.vector(mu_fout)) + 
    obs.clim[monstr.out]
  return(fi_ccr)
}
