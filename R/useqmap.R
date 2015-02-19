#' Quantile mapping using the qmap package
#' 
#' Computes bias correction with quantile mapping from the qmap package
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param ... additional arguments passed to \code{\link{fitQmap}} and
#'   \code{\link{doQmap}} (see \code{help(fitQmap)} and \code{help(doQmap)})
#'   
#' @details This function provides a wrapper to use the quantile mapping 
#'   functionality from the \code{qmap} package. All lead times are lumped 
#'   together, a potential seasonal cycle (lead-time dependency of quantile 
#'   correction) is not explicitly taken into account. Similarly, for the bias 
#'   corrected forecasts, all ensemble members are lumped together.
#'   
#' @examples
#' ## initialise forcast observation pairs
#' nens <- 51
#' signal <- outer(sin(seq(0,4,length=215)), sort(rnorm(30, sd=0.2)), '+')
#' fcst <- array(rgamma(length(signal)*nens, shape=3), c(dim(signal), nens)) +
#'   c(signal)
#' obs <- array(rnorm(length(signal), mean=2), dim(signal)) + 
#'   signal
#' fcst.debias <- biascorrection:::useqmap(fcst[,1:20,], 
#'   obs[,1:20], fcst.out=fcst[,21:30,], method='RQUANT', qstep=0.05, wet.day=FALSE)
#' oprob <- (seq(obs[,21:30]) - 1/3) / (length(obs[,21:30]) + 1/3)
#' plot(quantile(obs[,21:30], type=8, oprob), 
#'   quantile(fcst[,21:30,], type=8, oprob),
#'   type='l', lwd=2, xlab='Observed quantiles',
#'   ylab='Forecast quantiles',
#'   main='Out-of-sample validation for qqmap')
#' abline(c(0,1), lwd=2, lty=2)
#' lines(quantile(obs[,21:30], type=8, oprob),
#'   quantile(fcst.debias, type=8, oprob), lwd=2, col=2)
#' minprob <- min((11 - 1/3) / (length(obs[,1:20]) + 1/3), 0.05)
#' abline(v=quantile(obs[,1:20], type=8, prob=c(minprob, 1-minprob)), 
#'   lwd=2, lty=3)
#' text(quantile(obs[,1:20], type=8, 1-minprob) + 0.1, par('usr')[3] + 0.5, 
#'   'Extrapolated\ncorrection', adj=c(0,0), cex=0.67)
#' text(quantile(obs[,1:20], type=8, 1-minprob) - 0.1, par('usr')[3] + 0.5, 
#'   'Explicit quantile\ncorrection', adj=c(1,0), cex=0.67)
#' legend('topleft', c('No bias correction', 'useqmap (RQUANT)'), lwd=2, col=1:2, inset=0.05)
#' 
#' @keywords util
useqmap <- function(fcst, obs, fcst.out=fcst, ...){
  nens <- dim(fcst)[3]
  fq <- qmap::fitQmap(as.vector(obs), array(fcst, c(length(fcst)/nens, nens)), ...)
  fcst.debias <- array(qmap::doQmap(as.vector(fcst.out), fq, ...), dim(fcst.out))
  return(fcst.debias)
}