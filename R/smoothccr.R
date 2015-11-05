#' Daily CCR with monthly correction factors
#' 
#' Computes climate conserving recalibration on daily time series using
#' a linear regression approach with parameters following a polynomial (2nd order)
#' lead time dependency
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for lead time dependency
#' @param fcout.time forecast time for array to which bias correction is applied
#' for back compatibility with leave-one-out cross-validation in \code{\link{debias}}
#' @param span the parameter which controls the degree of smoothing (see \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction methods
#' 
#' @examples
#' signal <- rep(rnorm(30, sd=2), each=215)
#' fcst <- array(rnorm(30*215*51, mean=1, sd=rep(seq(0.5,2, length=30), each=215)), 
#' c(215, 30, 51)) + 0.5*sin(seq(0,4,length=215)) + signal*0.5
#' obs <- array(rnorm(30*215, mean=2), c(215, 30)) + sin(seq(0,4, length=215)) + signal
#' fc.time <- outer(1:215, 1981:2010, function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x)
#' fcst.debias <- biascorrection:::smoothccr(fcst, obs, fc.time=fc.time, span=0.5)
#' fcst.mon <- monmean(fcst, fc.time)
#' obs.mon <- monmean(obs, fc.time)
#' fcst.mondebias <- monmean(fcst.debias, fc.time)
#' 
#' @keywords util
smoothccr <- function(fcst, obs, fcst.out=fcst, fc.time, fcout.time=fc.time, 
                        span=min(1, 31/nrow(fcst)), ...){
  stopifnot(length(fcout.time) == length(fcst.out[,,1]))
  stopifnot(length(fc.time) >= length(fcst[,,1]))

  ## compute climatologies
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.out.ens <- rowMeans(fcst.out, dims=2, na.rm=T)

  in.day <- apply(fc.time, 2, function(x) x - x[1])
  out.day <- apply(fcout.time, 2, function(x) x - x[1])
  
  ## set up data frames
  in.df <- data.frame(fcst=c(fcst.ens), obs=c(obs), lead=c(in.day))
  out.df <- data.frame(fcst=c(fcst.out.ens), lead=c(out.day))
  
  ## fit regression
  flm <- lm(obs ~ fcst*poly(lead,2), in.df)
  
  ## predict the corrected forecast
  fpred <- array(predict(flm, newdata=out.df), dim(fcst.out.ens))
  
  ## add scaled anomalies
  fanom <- fcst.out - c(fcst.out.ens)
  res.sd <- sqrt(1 / ncol(obs) * apply(array(flm$res, dim(obs))**2, 1, sum))
  fout.sd <- sloess(sqrt(apply(apply(fanom, 1:2, sd, na.rm=T)**2, 1, mean)))
  fcst.debias <- fanom / fout.sd * res.sd + c(fpred)
  return(fcst.debias)  
}
