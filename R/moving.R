#'@name moving
#'  
#'@aliases moving_mul fastmoving fastmoving_mul
#'  
#'@title Moving Window Mean De-biasing
#'  
#'@description Computes moving window mean debiasing. A short-hand for multiplicative (\code{moving_mul}) 
#'  mean debiasing is also provided. Furthermore, a mean debiasing approach 
#'  in which the algorithm moves on in jumps rather than sequantially along lead
#'  times is provided via \code{fastmoving} and \code{fastmoving_mul}, 
#'  respectively.
#'  
#'@param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'  members
#'@param obs n x m matrix of veryfing observations
#'@param fcst.out array of forecast values to which bias correction should be 
#'  applied (defaults to \code{fcst})
#'@param window width of window of lead times to be used 
#'@param jump minimum number of lead times the moving window approach jumps
#'@param multiplicative logical, is mean debiasing to be applied 
#'  multiplicatively?
#'@param lower.bound is used to truncate output if set (e.g. to zero for 
#'  precipitation)
#'@param ... additional arguments for compatibility with other bias correction 
#'  methods
#'  
#'@details  
#'  The mean debiasing is lead time dependent, parameter \code{window} is used
#'  to select the number of lead times around the target lead time to be used for bias correction. For the begining
#'  and end of the series, the lead-time interval is kept constant, so that to
#'  estimate the correction for the first lead time, the first
#'  \code{window} lead times are used.
#'  
#' @examples
#' ## initialise forcast observation pairs
#' nens <- 51
#' signal <- outer(sin(seq(0,4,length=215)), sort(rnorm(30, sd=0.2)), '+') + 2
#' fcst <- list(raw=array(rgamma(length(signal)*nens, shape=2, scale=2), c(dim(signal), nens)) *
#'   c(signal) * (runif(length(signal)*nens) > 0.1))
#' obs <- array(rgamma(length(signal), shape=3, scale=1), dim(signal)) *
#'   signal * (runif(length(signal)) > 0.3)
#' fcst$moving <- biascorrection:::moving(fcst$raw[,1:20,], 
#'   obs[,1:20], fcst.out=fcst$raw, lower.bound=0)
#' fcst$moving_mul <- biascorrection:::moving_mul(fcst$raw[,1:20,], 
#'   obs[,1:20], fcst.out=fcst$raw, lower.bound=0)
#' oprob <- (seq(obs[,21:30]) - 1/3) / (length(obs[,21:30]) + 1/3)
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(density(obs[,21:30], from=0, to=80, bw=1), type='n',
#'      main='Distribution in validation period')
#' for (i in 1:length(fcst)) lines(density(fcst[[i]][,21:30,1], from=0, to=80, bw=1), lwd=2, lty=i)
#' lines(density(obs[,21:30], from=0, to=80, bw=1), lwd=2, col=2)
#' legend('topright', c('Observations', 'No bias correction', names(fcst)[-1]), 
#'        lwd=2, col=c(2,rep(1, length(fcst))), lty=c(1,seq(fcst)), inset=0.05)
#' plot(quantile(obs[,21:30], type=8, oprob), 
#'   quantile(fcst[[1]][,21:30,], type=8, oprob),
#'   type='l', lwd=2, xlab='Observed quantiles',
#'   ylab='Forecast quantiles',
#'   main='Out-of-sample validation for moving')
#' abline(c(0,1), lwd=2, col=2)
#' for (i in 2:length(fcst)) lines(quantile(obs[,21:30], type=8, oprob),
#'   quantile(fcst[[i]][,21:30,], type=8, oprob), lwd=2, lty=i)
#' legend('topleft', c('No bias correction', names(fcst)[-1]), lwd=2, lty=seq(along=fcst), inset=0.05)
#' par(oldpar)
#' 
#'@keywords util
moving <- function(fcst, obs, fcst.out=fcst, 
                  window=min(nrow(fcst), 31), jump=1,
                  multiplicative=FALSE, lower.bound=NULL, 
                  ...){

  ## check input data
  stopifnot(is.matrix(obs), is.array(fcst), dim(fcst)[1:2] == dim(obs))
  ## number of lead times
  nlead <- nrow(fcst)
  ## set everything to missing if fcst / obs is missing
  fcst[rep(is.na(obs), length=length(fcst))] <- NA
  ## multiply obs to get corresponding array (in case of missing values)
  obs <- array(obs, dim(fcst))
  obs[is.na(fcst)] <- NA
  ## compute 'ensemble mean'
  fcst.mn <- rowMeans(fcst, dims=2, na.rm=T)
  obs.mn <- rowMeans(obs, dims=2, na.rm=T)

  ## check window and jump
  stopifnot(window > 0, jump <= window)
  window <- min(window, nlead)
  

  ## start main de-biasing algorithm

  ## initialize output
  fcst.debias <- NA*fcst.out
  
  ## figure out number of jumps
  njump <- ceiling((nlead - window) / jump) + 1
  jrest <- ceiling((window - jump) / 2)
  ## new function to simplify notation
  jseq <- function(from=1) seq(from=from, by=jump, length=njump)
  mini <- pmin(jseq(1), nlead - window + 1)
  maxi <- pmin(jseq(window), nlead)
  mino <- pmin(jseq(1 + jrest), nlead - window + jrest + 1)
  mino[1] <- 1
  maxo <- pmin(jseq(jrest + jump), nlead)
  maxo[njump] <- nlead
  
  ## loop through lead times
  for (i in seq(along=mini)){
    ## lead times for estimating correction
    ind <- seq(mini[i], maxi[i])
    ## lead times to be bias corrected
    ind2 <- seq(mino[i], maxo[i])
    
    ## compute bias
    f.mn <- mean(fcst.mn[ind,], na.rm=T)
    o.mn <- mean(obs.mn[ind,], na.rm=T)
    if (multiplicative){
      fcorr <- ifelse(f.mn == 0, 1, o.mn / f.mn)
      fcst.debias[ind2,,] <- fcst.out[ind2,,]*fcorr  
    } else {
      fcst.debias[ind2,,] <- fcst.out[ind2,,] - f.mn + o.mn
    }
  } ## end of loop on lead times
  

  if (!is.null(lower.bound)){
    fcst.debias[fcst.debias < lower.bound] <- lower.bound
  }
  
  return(fcst.debias)
}


#'@rdname moving
#'
moving_mul <- function(..., multiplicative=TRUE){
  moving(..., multiplicative=multiplicative)
}

#'@rdname moving
#'
fastmoving <- function(fcst, window=min(nrow(fcst), 31), jump=min(window, 11), ...){
  moving(fcst=fcst, window=window, jump=jump, ...)
}

#'@rdname moving
#'
fastmoving_mul <- function(fcst, window=min(nrow(fcst), 31), 
                          jump=min(window, 11), multiplicative=TRUE, ...){
  fastmoving(fcst=fcst, window=window, jump=jump, multiplicative=multiplicative, ...)
}
