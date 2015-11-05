#' Daily CCR with monthly correction factors
#' 
#' Computes climate conserving recalibration on smoothed monthly means
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for monthly aggregation
#' @param fcout.time forecast time for array to which bias correction is applied
#' for back compatibility with leave-one-out cross-validation in \code{\link{debias}}
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @examples
#' fcst <- array(rnorm(30*215*51, mean=1, sd=rep(seq(0.5,2, length=30), each=215)), 
#' c(215, 30, 51)) + 0.5*sin(seq(0,4,length=215))
#' obs <- array(rnorm(30*215, mean=2), c(215, 30)) + sin(seq(0,4, length=215))
#' fc.time <- outer(1:215, 1981:2010, function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x)
#' fcst.debias <- biascorrection:::ccr_monthly(fcst, obs, fc.time=fc.time, span=0.5)
#' fcst.mon <- monmean(fcst, fc.time)
#' obs.mon <- monmean(obs, fc.time)
#' fcst.mondebias <- monmean(fcst.debias, fc.time)
#' 
#' @keywords util
ccr_monthly <- function(fcst, obs, fcst.out=fcst, fc.time, fcout.time=fc.time, 
                        span=min(1, 31/nrow(fcst)), ...){
  stopifnot(length(fcout.time) == length(fcst.out[,,1]))
  stopifnot(length(fc.time) >= length(fcst[,,1]))
  
  ## compute monthly forecasts and observations
  fcst.mon <- monmean(fcst, fc.time[seq(along=fcst[,,1])])
  obs.mon <- monmean(obs, fc.time[seq(obs)])
  
  ## compute climatologies
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  fcst.clim <- sloess(fcst.mn, span=span)
  obs.clim <- sloess(obs.mn, span=span)  
  fcst.monclim <- rowMeans(fcst.mon, dims=1, na.rm=T)
  obs.monclim <- rowMeans(obs.mon, dims=1, na.rm=T)
  ## additional monthly bias to assure monthly means fit
  monbias <- rowMeans(monmean(fcst.ens - obs - (fcst.clim - obs.clim), 
                            fc.time[seq(along=fcst.ens)]), 
                      dims=1, na.rm=T)
  
  ## compute ccr_monthly prerequisites (Notation as in Weigel et al. 2009)
  x <- obs.mon - obs.monclim
  sig_x <- sqrt(apply(x**2, 1, mean))
  fi <- fcst.mon - fcst.monclim
  mu_f <- rowMeans(fi, dims=2)
  sig_mu <- sqrt(apply(mu_f**2, 1, mean))
  sig_ens <- sqrt(apply(apply(fi, 1:2, sd)**2, 1, mean))
  rho <- diag(cor(t(mu_f), t(x)))
  rr <- rho*sig_x / sig_mu
  ss <- sqrt(1 - rho**2) * sig_x / sig_ens
  
  ## put everything back together
  monstr.out <- format(fcout.time, '%m')
  fi_out <- fcst.out - fcst.clim
  mu_fout <- rowMeans(fi_out, dims=2, na.rm=T)
  mu_fout.mon <- monmean(mu_fout, fcout.time)
  ## fix to be able to get the correct index out
  colnames(mu_fout.mon) <- as.character(1:ncol(mu_fout.mon))
  ind <- cbind(monstr.out, rep(1:ncol(mu_fout.mon), each=nrow(mu_fout)))
  fi_ccr <- (rr*mu_fout.mon)[ind] + 
    ss[monstr.out]*(fi_out - as.vector(mu_fout)) + 
    obs.clim - monbias[monstr.out]
  return(fi_ccr)
}
