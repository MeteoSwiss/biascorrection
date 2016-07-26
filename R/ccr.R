#' Climate Conserving Recalibration
#' 
#' Climate conserving recalibration with correction factors for each lead time 
#' without smoothing (application to daily series not recommended).
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param type if set to "prediction", additional inflation of the spread as in 
#'   linear regression (see details)
#' @param ... additional arguments for compatibility with other bias correction 
#'   methods
#'   
#' @details In \code{calibration} mode, the climate conserving recalibration 
#'   (CCR) follows Weigel \emph{et al.} (2008). In \code{prediction} mode, the 
#'   CCR spread correction is expanded to take into account additional 
#'   uncertainty from the signal calibration following linear regression theory.
#'   I.e. the inflation factor \code{s.pred} is
#'   
#'   \deqn{s.pred = s \sqrt{1 + 1/n + f_0^2 / \sum{f_j^2}}}{s.pred = s * sqrt(1
#'   + 1/n + fo^2 / sum(fj^2))}
#'   
#'   where \eqn{f_0}{fo} is the ensemble mean forecast anomaly that is to be 
#'   adjusted, and \eqn{\sum{f_j^2}}{sum(fj^2)} is the sum of the squared 
#'   ensemble mean forecast anomalies in the calibration set. \eqn{n} is the 
#'   number of forecast instances in the calibration set.
#'   
#' @references Weigel, A., M. Liniger and C. Appenzeller (2008). Seasonal 
#'   Ensemble Forecasts: Are Recalibrated Single Models Better than Multimodels?
#'   \emph{Monthly Weather Review}, 137(4), 1460-1479.
#'   
#' @examples
#' fcst <- array(rnorm(3000*1*51, mean=1, sd=rep(seq(0.5,2, length=3000), each=1)), 
#' c(1, 3000, 51)) + 0.5*sin(seq(0,4,length=1))
#' obs <- array(rnorm(3000, mean=2), c(1, 3000)) + sin(seq(0,4, length=1))
#' fcst.debias <- biascorrection:::ccr(fcst, obs)
#' f.rmse <- sqrt(apply((obs - apply(fcst.debias, 2, mean))**2, 1, mean))
#' f.sd <- sqrt(mean(apply(fcst.debias, 2, sd)**2))
#' f.sd / f.rmse ## should be exactly 1
#' mean(fcst.debias - obs[1,]) ## should be 0 (rounding errors)
#' 
#' @keywords util
ccr <- function(fcst, obs, fcst.out=fcst, type=c("calibration", "prediction"), ...){
  type <- match.arg(type)
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  ## compute climatology
  fcst.clim <- rowMeans(fcst.ens, dims=1, na.rm=T)
  obs.clim <- rowMeans(obs, dims=1, na.rm=T)
  
  ## compute CCR prerequisites (Notation as in Weigel et al. 2009)
  x <- obs - obs.clim
  sig_x <- sqrt(apply(x**2, 1, mean))
  fi <- fcst - fcst.clim
  mu_f <- rowMeans(fi, dims=2)
  sig_mu <- sqrt(apply(mu_f**2, 1, mean))
  sig_ens <- sqrt(apply(apply(fi, 1:2, sd)**2, 1, mean))
  rho <- diag(cor(t(mu_f), t(x)))
  rho[is.na(rho)] <- 0
  rr <- rho*sig_x / sig_mu
  rr[sig_mu == 0] <- 0
  ss <- sqrt(1 - rho**2) * sig_x / pmax(sig_ens, 1e-12)
  ## put everything back together (with de-biasing)
  fi_out <- fcst.out - fcst.clim
  mu_fout <- rowMeans(fi_out, dims=2, na.rm=T)
  ## add additional spread correction for out-of-sample calibration  
  if (type == 'prediction'){
    ss <- ss*sqrt(1 + 1/ncol(mu_f) + mu_fout**2 / apply(mu_f**2, 1, sum))
  }
  fi_ccr <- as.vector(rr * mu_fout) + c(ss) * (fi_out - as.vector(mu_fout)) + obs.clim  
  return(fi_ccr)
}
