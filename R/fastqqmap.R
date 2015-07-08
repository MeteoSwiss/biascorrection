#'@name
#'fastqqmap
#'
#'@aliases
#'fastqqmap_mul
#'
#'@title
#'Quantile mapping
#'
#'@description
#'Computes bias correction with quantile mapping
#'
#'@param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'  members
#'@param obs n x m matrix of veryfing observations
#'@param fcst.out array of forecast values to which bias correction should be 
#'  applied (defaults to \code{fcst})
#'@param multiplicative logical, is quantile correction to be applied 
#'  multiplicatively?
#'@param lower.bound is used to truncate output if set (e.g. to zero for
#'  precipitation)
#'@param anomalies logical, should quantile mapping be applied to forecast and 
#'  observed anomalies (from forecast ensemble mean) only?
#'@param window width of window to be used for quantile mapping. To increase 
#'  computation speed, the window moves along at approximately 1/6 of the window
#'  width (i.e. jumps by 15 days for a 91 day window)
#'@param ... additional arguments for compatibility with other bias correction 
#'  methods
#'  
#'@details The quantile mapping algorithm estimates quantile correction factors 
#'  for \code{n} quantiles. For each forecast value in \code{fcst.out}, the 
#'  percentile within which the value falls in the distribution of input 
#'  forecasts \code{fcst} is determined and the corresponding quantile 
#'  correction applied. For multiplicative quantile mapping 
#'  (\code{multiplicative = TRUE}), the bias corrected forecast 
#'  (\code{fcst.out}) is divided by the ratio of forecast to observed quantiles,
#'  whereas for additive quantile mapping \code{multiplicative = FALSE}, the 
#'  difference between the forecast and observed quantiles are subtracted from 
#'  \code{fcst.out}. The quantiles are estimated for at least 50 discrete values
#'  from the 5th to 95th percentile, or, if there are enough observations for 
#'  \code{n = n_obs / 10} discrete quantiles excluding the 5 smallest and 
#'  largest values. If \code{anomalies} is set, forecast and observed anomalies 
#'  are computed with reference to the forecast ensemble mean (the signal) and 
#'  the quantile mapping is only applied to the anomalies with the signal being 
#'  left uncorrected.
#'  
#'  The quantile mapping is lead time dependent, parameter \code{windo} is used
#'  to select the number of lead times to be included in the quantile
#'  estimation. For the beginning and end of the series, the lead-time interval
#'  is kept constant, so that to estimate the quantile correction for the first
#'  lead time, the first \code{window} lead times are used.
#'  
#'  
#' @examples
#' ## initialise forcast observation pairs
#' nens <- 51
#' signal <- outer(sin(seq(0,4,length=215)), sort(rnorm(30, sd=0.2)), '+') + 2
#' fcst <- array(rgamma(length(signal)*nens, shape=2, scale=2), c(dim(signal), nens)) *
#'   c(signal) * (runif(length(signal)*nens) > 0.1)
#' obs <- array(rgamma(length(signal), shape=3, scale=1), dim(signal)) *
#'   signal * (runif(length(signal)) > 0.3)
#' fcst.debias <- biascorrection:::fastqqmap(fcst[,1:20,], 
#'   obs[,1:20], fcst.out=fcst[,21:30,], lower.bound=0, multiplicative=TRUE)
#' oprob <- (seq(obs[,21:30]) - 1/3) / (length(obs[,21:30]) + 1/3)
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(density(fcst.debias, from=0, to=80), lwd=2, col=1, 
#'      main='Distribution in validation period')
#' lines(density(fcst[,21:30,], from=0, to=80), lwd=2, lty=2)
#' lines(density(obs[,21:30], from=0, to=80), lwd=2, col=2)
#' legend('topright', c('Observations', 'No bias correction', 'fastqqmap'), 
#'        lwd=2, col=c(2,1,1), lty=c(1,2,1), inset=0.05)
#' plot(quantile(obs[,21:30], type=8, oprob), 
#'   quantile(fcst[,21:30,], type=8, oprob),
#'   type='l', lwd=2, lty=2, xlab='Observed quantiles',
#'   ylab='Forecast quantiles',
#'   main='Out-of-sample validation for fastqqmap')
#' abline(c(0,1), lwd=2, col=2)
#' lines(quantile(obs[,21:30], type=8, oprob),
#'   quantile(fcst.debias, type=8, oprob), lwd=2)
#' legend('topleft', c('No bias correction', 'fastqqmap'), lwd=2, lty=2:1, inset=0.05)
#' par(oldpar)
#' 
#'@keywords util
fastqqmap <- function(fcst, obs, fcst.out=fcst, anomalies=FALSE, multiplicative=FALSE, lower.bound=NULL, window=min(nrow(fcst), 91), ...){
  ## estimate the quantile correction excluding the 5 smallest/largest values
  ## minprob <- min((6 - 1/3) / (ncol(obs)* window + 1/3), 0.05)
  ## estimate the quantile correction for the full range
  minprob <- min((2/3) / (ncol(obs)* window + 1/3), 0.05)
  prob <- seq(minprob, 1 - minprob, length=floor(max(50, ncol(obs)*window/20)))
  ## fq <- quantile(fcst, type=8, prob=prob)
  if (anomalies){
    fcst.ens <- rowMeans(fcst, dims=2)
    fcst.out.ens <- rowMeans(fcst.out, dims=2)
    if (multiplicative){
      fcst.anom <- fcst / c(fcst.ens)
      obs.anom <- obs/fcst.ens
      fcst.out.anom <- fcst.out / c(fcst.out.ens)
    } else {
      fcst.anom <- fcst - c(fcst.ens)
      obs.anom <- obs - fcst.ens
      fcst.out.anom <- fcst.out - c(fcst.out.ens)
    }    
  } else {
    fcst.anom <- fcst
    obs.anom <- obs
    fcst.out.anom <- fcst.out
  }
  nlead <- nrow(fcst)
  if (window >= nlead){
    fq <- rowMeans(apply(fcst.anom, 3, quantile, type=8, prob=prob))
    oq <- quantile(obs.anom, type=8, prob=prob)
    ## find boundaries in between quantiles
    fqbnds <- fq[-length(fq)] + 0.5*diff(fq)
    fout.qi <- findInterval(fcst.out.anom, fqbnds) + 1
    if (multiplicative){
      fcst.debias <- fcst.out * (oq / fq)[fout.qi]
    } else {
      fcst.debias <- fcst.out + (oq - fq)[fout.qi]
    }
  } else {
    ## do everything along the lead times
    fcst.debias <- NA*fcst.out.anom
    ## figure out jump width
    njump <- floor(window/6)
    isteps <- seq(ceiling(window/2),
                  nlead - floor(window/2),
                  by=njump)
    for (i in seq(along=isteps)){
      ## lead times for estimating quantiles
      ind <- seq(if(i == 1) 1 else isteps[i] - floor(window/2), 
                 if (i == length(isteps)) nlead else isteps[i] + floor(window/2))
      ## lead times to be bias corrected
      ind2 <- seq(if(i == 1) 1 else isteps[i] - floor(njump/2), 
                  if (i == length(isteps)) nlead else isteps[i] + floor(njump/2))
      oq <- quantile(obs.anom[ind,], prob=prob, type=8)
      fq <- quantile(fcst.anom[ind,,], prob=prob, type=8)
      ## fq <- rowMeans(apply(fcst.anom[ind,,], 3, quantile, prob=prob, type=8))
      ## find boundaries in between quantiles
      fqbnds <- fq[-length(fq)] + 0.5*diff(fq)
      
      ## get the probability of the output fcst given the forecast
      ## i.e. the reverse quantile function
      ## assume constant correction by discrete quantiles
      fout.qi <- findInterval(fcst.out.anom[ind2,,], fqbnds) + 1
      if (multiplicative){
        qcorr <- oq/fq
        qcorr[fq == 0] <- 1
        fcst.debias[ind2,,] <- fcst.out[ind2,,] * qcorr[fout.qi]
      } else {
        fcst.debias[ind2,,] <- fcst.out[ind2,,] - (fq - oq)[fout.qi]    
      }      
    } ## end of loop on lead times
  } ## end of if on window length
  if (!is.null(lower.bound)){
    fcst.debias[fcst.debias < lower.bound] <- lower.bound
  }
  return(fcst.debias)
}

#'@rdname fastqqmap
#'@export
#'
fastqqmap_mul <- function(fcst, obs, fcst.out=fcst, anomalies=FALSE, multiplicative=TRUE, lower.bound=NULL, window=min(nrow(fcst), 91), ...){
  fastqqmap(fcst=fcst, obs=obs, fcst.out=fcst.out, 
            anomalies=anomalies, multiplicative=multiplicative,
            lower.bound=lower.bound, window=window, ...=...)
}
