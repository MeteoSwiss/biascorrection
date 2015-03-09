#' Quantile mapping using the qmap package
#' 
#' Computes bias correction with quantile mapping from the qmap package
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param nn number of lead times to be included around the lead time to be 
#'   calibrated (i.e. \code{ceiling((nn - 1)/2)} before and \code{floor((nn - 
#'   1)/2)} after). The interval to estimate quantile corrections is held 
#'   constant, that is, the first \code{nn} lead times are used to estimate the 
#'   quantile correction for lead time one. The default is to lump all lead
#'   times together \code{nn = nrow(fcst)}.
#' @param exact logical, should surrounding lead-times be used for quantiles of 
#'   forecast as well?
#' @param ... additional arguments passed to \code{\link{fitQmap}} and 
#'   \code{\link{doQmap}} (see \code{help(fitQmap)} and \code{help(doQmap)})
#'   
#' @details This function provides a wrapper to use the quantile mapping 
#'   functionality from the \code{qmap} package. All lead times are lumped 
#'   together, a potential seasonal cycle (lead-time dependency of quantile 
#'   correction) is not explicitly taken into account. Similarly, for the bias 
#'   corrected forecasts, all ensemble members are lumped together.
#'   
#'   The quantile mapping is lead time dependent, parameter \code{nn} is used to
#'   select the number of lead times on either side of the lead time that is to 
#'   be corrected to be included in the quantile estimation. For the beginning 
#'   and end of the series, the lead-time interval is kept constant, so that to 
#'   estimate the quantile correction for the first lead time, the first 
#'   \code{nn} lead times are used. If \code{exact = FALSE}, the lead time 
#'   dependent quantiles for the forecast are directly estimated from single
#'   lead times without the surrounding \code{2*nn} lead times. Note, however,
#'   that unlike for \code{\link{qqmap}}, this does not speed up computation.
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
#' minprob <- 0.05
#' abline(v=quantile(obs[,1:20], type=8, prob=c(minprob, 1-minprob)), 
#'   lwd=2, lty=3)
#' text(quantile(obs[,1:20], type=8, 1-minprob) + 0.1, par('usr')[3] + 0.5, 
#'   'Extrapolated\ncorrection', adj=c(0,0), cex=0.67)
#' text(quantile(obs[,1:20], type=8, 1-minprob) - 0.1, par('usr')[3] + 0.5, 
#'   'Explicit quantile\ncorrection', adj=c(1,0), cex=0.67)
#' legend('topleft', c('No bias correction', 'useqmap (RQUANT)'), 
#'   lwd=2, col=1:2, inset=0.05)
#' 
#' @keywords util
useqmap <- function(fcst, obs, fcst.out=fcst, nn=nrow(fcst), exact=FALSE, ...){
  nens <- dim(fcst)[3]
  nlead <- nrow(fcst)
  if (nn >= nlead){
    fq <- qmap::fitQmap(as.vector(obs), array(fcst, c(length(fcst)/nens, nens)), ...)
    fcst.debias <- array(qmap::doQmap(as.vector(fcst.out), fq, ...), dim(fcst.out))    
  } else if (exact){
    fcst.debias <- fcst.out*NA
    indold <- rep(1, nn)
    for (i in 1:nlead){
      ind <- seq(0, nn - 1) + min(max(i - ceiling((nn - 1)/2), 1), nlead - nn + 1)
      if (! all(ind == indold)){
        fq <- qmap::fitQmap(as.vector(obs[ind,]), matrix(fcst[ind,,], nn*ncol(fcst), nens), ...)
      } 
      indold <- ind
      fcst.debias[i,,] <- c(qmap::doQmap(as.vector(fcst.out[i,,]), fq, ...))
    }
  } else {
    fcst.debias <- fcst.out*NA
    indold <- rep(1, nn)
    for (i in 1:nlead){
      ind <- seq(0, nn - 1) + min(max(i - ceiling((nn - 1)/2), 1), nlead - nn + 1)
      fq <- qmap::fitQmap(as.vector(obs[ind,]), fcst[i,,], ...)
      fcst.debias[i,,] <- c(qmap::doQmap(as.vector(fcst.out[i,,]), fq, ...))
    }
  }
  return(fcst.debias)
}
