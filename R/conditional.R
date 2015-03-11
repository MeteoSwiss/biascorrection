#' Bias conditional on signal
#' 
#' Compute calibration for biases that are conditional on signal
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param span the parameter which controls the degree of smoothing (see 
#'   \code{\link{loess}})
#' @param ... additional arguments for compatibility with other bias correction 
#'   methods
#'   
#' @details This bias correction method assumes that the bias can be decomposed 
#'   into a stationary component and a component linearly related to the 
#'   ensemble- average prediction (i.e. the magnitude of the bias for forecasts 
#'   that are farther from the mean). Let \eqn{X_{i,k}}{Xik} be forecast member 
#'   \eqn{k} at forecast instance \eqn{i}, and \eqn{\bar{X_i}}{mn.Xi} the 
#'   ensemble mean for instance \eqn{i}, and \eqn{\bar{X}}{mn.X} the mean 
#'   forecast averaged over all forecast instances. Given verifying observations
#'   \eqn{y_i}{yi}, the bias corrected forecast \eqn{Z_{i,k}}{Zik} is
#'   
#'   \deqn{Z_{i,k} = X_{i,k} - b{i}}{Zik = Xik - bi}
#'   
#'   with the bias \eqn{b_{i}}{bi} estimated from the forecast anomalies 
#'   \eqn{X'_i = \bar{X_i} - \bar{X}}{X'i = mn.Xi - mn.X} and observed anomalies
#'   \eqn{y'_i = y_i - \bar{y}}{y'i = yi - mn.y}
#'   
#'   \deqn{b_{i} = E[X'_i - y'_i | X'_i] + \bar{X} - \bar{y}}{bi = E[X'i - y'i |
#'   X'i] + mn.X - mn.y}
#'   
#'   In the above, the index for different forecast lead times are dropped. In
#'   practice, the mean forecast \eqn{\bar{X}}{mn.X} and observations
#'   \eqn{\bar{y}}{mn.y} may be vectors with entries for different lead times.
#'   The averaged quantities are smoothed to result in a smooth climatology for
#'   different lead times.
#'   
#' @examples
#' ## initialise forcast observation pairs
#' seasonal <- sin(seq(0,4,length=215))
#' signal <- outer(seasonal, rnorm(30), '+')
#' fcst <- array(rnorm(215*30*51), c(215, 30, 15)) + 
#'   2*c(signal)
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) +
#'   signal
#' fcst.debias <- biascorrection:::conditional(fcst[,1:20,], 
#'   obs[,1:20], fcst.out=fcst[,21:30,], span=0.5)
#' fcst.debias2 <- biascorrection:::smooth(fcst[,1:20,], 
#'   obs[,1:20], fcst.out=fcst[,21:30,], span=0.5)
#' boxplot(cbind(raw=c(fcst[,21:30,]) - c(obs[,21:30]),
#'   smooth=c(fcst.debias2) - c(obs[,21:30]),
#'   conditional=c(fcst.debias) - c(obs[,21:30])),
#'   ylab='distribution of forecast errors',
#'   main='Out-of-sample validation for conditional bias')
#' abline(h=0, lwd=2, lty=2)
#' 
#' @keywords util
conditional <- function(fcst, obs, fcst.out=fcst, span=min(1, 31/nrow(fcst)), ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  fcst.clim <- sloess(fcst.mn, span=span)
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- sloess(obs.mn, span=span)
  
  fcst.out.ens <- rowMeans(fcst.out, dims=2)
  
  in.df <- data.frame(fcst=c(fcst.ens - fcst.clim),
                      obs=c(obs - obs.clim))
  out.df <- data.frame(fcst=c(fcst.out.ens - fcst.clim))
  
  f.lm <- lm(obs ~ fcst - 1, in.df)
  fcst.debias <- fcst.out - c(fcst.out.ens) + 
    predict(f.lm, newdata=out.df) + obs.clim
  
  ## old approach 
  ## residual anomalies independent of climatology
  ## anom <- fcst.ens - fcst.clim - (obs - obs.clim)
  
  ## compute conditional bias (subtracting constant term)
  ## anom.lm <- lm(anom ~ fcst - 1, data.frame(anom=c(anom), 
  ##                                           fcst=c(fcst.ens - fcst.clim)))

  ## compute debiased forecast
  ## fcst.out.ens <- rowMeans(fcst.out, dims=2)
  ## fcst.debias <- fcst.out - 
  ##   predict(anom.lm, newdata=data.frame(fcst=c(fcst.out.ens - fcst.clim))) +
  ##   obs.clim - fcst.clim
  
  return(fcst.debias)
}

