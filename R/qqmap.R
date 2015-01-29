#' Quantile mapping
#' 
#' Computes bias correction with quantile mapping
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param multiplicative logical, is quantile correction to be applied 
#'   multiplicatively?
#' @param ... additional arguments for compatibility with other bias correction 
#'   methods
#'   
#' @details The quantile mapping algorithm estimates quantile correction factors
#'   for \code{n} quantiles. For each forecast value in \code{fcst.out}, the 
#'   percentile within which the value falls in the distribution of input 
#'   forecasts \code{fcst} is determined and the corresponding quanile 
#'   correction applied. For multiplicative quantile mapping 
#'   (\code{multiplicative = TRUE}), the bias corrected forecast 
#'   (\code{fcst.out}) is divided by the ratio of forecast to observed 
#'   quantiles, whereas for additive quantile mapping \code{multiplicative = 
#'   FALSE}, the difference between the forecast and observed quantiles are 
#'   subtracted from \code{fcst.out}. The quantiles are estimated for at least
#'   100 discrete values from the 5th to 95th percentile, or, if there are
#'   enough observations for \code{n = n_obs / 10} discrete quantiles excluding
#'   the 10 smallest and largest values.
#'   
#' @note The quantile mapping provided here does not take into account lead-time
#'   dependent quantile-quantile relationships. Instead, all lead times are 
#'   lumped together.
#'   
#' @examples
#' ## initialise forcast observation pairs
#' nens <- 51
#' signal <- outer(sin(seq(0,4,length=215)), sort(rnorm(30, sd=0.2)), '+')
#' fcst <- array(rgamma(length(signal)*nens, shape=3), c(dim(signal), nens)) +
#'   c(signal)
#' obs <- array(rnorm(length(signal), mean=2), dim(signal)) + 
#'   signal
#' fcst.debias <- biascorrection:::qqmap(fcst[,1:20,], 
#'   obs[,1:20], fcst.out=fcst[,21:30,])
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
#' legend('topleft', c('No bias correction', 'qqmap'), lwd=2, col=1:2, inset=0.05)
#' 
#' @keywords util
qqmap <- function(fcst, obs, fcst.out=fcst, multiplicative=FALSE, ...){
  ## only estimate the quantile correction from 5 to 95th percentile
  ## or excluding the 10 smallest and largest values
  minprob <- min((11 - 1/3) / (length(obs) + 1/3), 0.05)
  prob <- seq(minprob, 1 - minprob, length=max(100, length(obs)/10))
  ## fq <- quantile(fcst, type=8, prob=prob)
  fq <- rowMeans(apply(fcst, 3, quantile, type=8, prob=prob))
  oq <- quantile(obs, type=8, prob=prob)
  
  ## find boundaries in between quantiles
  fqbnds <- fq[-length(fq)] + 0.5*diff(fq)
  
  ## get the probability of the output fcst given the forecast
  ## i.e. the reverse quantile function
  ## assume constant correction by discrete quantiles
  fout.qi <- findInterval(fcst.out, fqbnds) + 1
  if (multiplicative){
    fcst.debias <- fcst.out / (fq/oq)[fout.qi]
  } else {
    fcst.debias <- fcst.out - (fq - oq)[fout.qi]    
  }
  
  return(fcst.debias)
}
