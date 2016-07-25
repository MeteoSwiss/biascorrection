#'@name qqmap
#'  
#'@aliases qqmap_mul fastqqmap fastqqmap_mul
#'  
#'@title Quantile Mapping
#'  
#'@description Computes bias correction with quantile mapping (i.e. additive 
#'  quantile correction). A short-hand for multiplicative (\code{qqmap_mul}) 
#'  quantile mapping is also provided. Furthermore, a quantile mapping approach 
#'  in which the algorithm moves on in jumps rather than sequantially along lead
#'  times is provided via \code{fastqqmap} and \code{fastqqmap_mul}, 
#'  respectively.
#'  
#'@param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'  members
#'@param obs n x m matrix of veryfing observations
#'@param fcst.out array of forecast values to which bias correction should be 
#'  applied (defaults to \code{fcst})
#'@param minprob probability boundary for quantiles. Quantiles from 
#'  \code{minprob} to \code{1 - minprob} are estimated. Defaults to \code{0.01} 
#'  for percentiles from 1 to 99.
#'@param nprob number of quantiles to estimate correction factors. By default 
#'  \code{min((1 - minprob) / minprob, 201)} quantiles are estimated.
#'@param window width of window to be used for quantile mapping
#'@param jump minimum number of days the moving quantile window jumps (see 
#'  below)
#'@param multiplicative logical, is quantile correction to be applied 
#'  multiplicatively?
#'@param lower.bound is used to truncate output if set (e.g. to zero for 
#'  precipitation)
#'@param debias logical, should quantile mapping be applied to anomalies from 
#'  smoothed climatology (only additive)?
#'@param smoothobs logical, should observation climatology be smoothed?
#'@param smooth logical, should forecast climatology be smoothed?
#'@param span the smoothing bandwidth (see \code{\link{loess}})
#'@param anomalies logical, should quantile mapping be applied to forecast and 
#'  observed anomalies (from forecast ensemble mean) only (see below)?
#'@param ... additional arguments for compatibility with other bias correction 
#'  methods
#'  
#'@details The quantile mapping algorithm estimates quantile correction factors 
#'  for \code{q} quantiles. For each forecast value in \code{fcst.out}, the 
#'  percentile within which the value falls in the distribution of input 
#'  forecasts \code{fcst} is determined and the corresponding quantile 
#'  correction applied. For multiplicative quantile mapping 
#'  (\code{multiplicative = TRUE}), the bias corrected forecast 
#'  (\code{fcst.out}) is divided by the ratio of forecast to observed quantiles,
#'  whereas for additive quantile mapping \code{multiplicative = FALSE} (the
#'  default), the difference between the forecast and observed quantiles are
#'  subtracted from \code{fcst.out}.
#'  
#'  The quantile mapping is lead time dependent, parameter \code{window} is used
#'  to select the number of lead times on either side of the lead time that is 
#'  to be corrected to be included in the quantile estimation. For the begining
#'  and end of the series, the lead-time interval is kept constant, so that to
#'  estimate the quantile correction for the first lead time, the first
#'  \code{window} lead times are used. If \code{exact = FALSE}, the lead time
#'  dependent quantiles for the forecast are directly estimated from single lead
#'  times without the surrounding \code{window} lead times. This is a quick and
#'  dirty fix to speed up processing.
#'  
#'@section De-biasing: If \code{debias = TRUE}, the quantile correction is 
#'  applied to the anomalies from the long-term lead-time dependent climatology 
#'  (from the observations and forecasts respectively). The quantile corrected 
#'  forecast anomalies are finally added to the observed climatology to produce 
#'  an approximately unbiased, quantile corrected forecast. If \code{smoothobs =
#'  TRUE} and/or \code{smooth = TRUE}, the lead-time dependent climatology of 
#'  the observations and/or forecasts are smoothed using a 
#'  \code{\link[stats]{loess}} smoother with bandwidth \code{span}. Whether 
#'  there are use-cases for which such an \emph{a priori} de-biasing is 
#'  beneficial is not obvious and needs further exploration.
#'  
#'@section Anomalies: If \code{anomalies} is set, forecast and observed 
#'  anomalies are computed with reference to the forecast ensemble mean (the 
#'  signal) and the quantile correction is only applied to the anomalies with 
#'  the signal being left uncorrected. It is speculated that such an approach 
#'  may more explicitly reflect the skill / calibration relationship of 
#'  forecasts (i.e. the option has been implemented without theoretical 
#'  underpinning nor an existing use-case that illustrates the benefits of 
#'  applying quantile corrections as described in this paragraph ;-).
#'  
#' @examples
#' ## initialise forcast observation pairs
#' nens <- 51
#' signal <- outer(sin(seq(0,4,length=215)), sort(rnorm(30, sd=0.2)), '+') + 2
#' fcst <- list(raw=array(rgamma(length(signal)*nens, shape=2, scale=2), c(dim(signal), nens)) *
#'   c(signal) * (runif(length(signal)*nens) > 0.1))
#' obs <- array(rgamma(length(signal), shape=3, scale=1), dim(signal)) *
#'   signal * (runif(length(signal)) > 0.3)
#' fcst$qqmap <- biascorrection:::qqmap(fcst$raw[,1:20,], 
#'   obs[,1:20], fcst.out=fcst$raw, lower.bound=0)
#' fcst$qqmap_mul <- biascorrection:::qqmap_mul(fcst$raw[,1:20,], 
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
#'   main='Out-of-sample validation for qqmap')
#' abline(c(0,1), lwd=2, col=2)
#' for (i in 2:length(fcst)) lines(quantile(obs[,21:30], type=8, oprob),
#'   quantile(fcst[[i]][,21:30,], type=8, oprob), lwd=2, lty=i)
#' legend('topleft', c('No bias correction', names(fcst)[-1]), lwd=2, lty=seq(along=fcst), inset=0.05)
#' par(oldpar)
#' 
#'@keywords util
qqmap <- function(fcst, obs, fcst.out=fcst, 
                  minprob=0.01, nprob=min((1 - minprob)/minprob, 201),
                  window=min(nrow(fcst), 31), jump=1,
                  multiplicative=FALSE, lower.bound=NULL, 
                  anomalies=FALSE, debias=FALSE,
                  smoothobs=TRUE, smooth=smoothobs,
                  span=min(91/nrow(fcst), 1), ...){

  ## check input data
  stopifnot(is.matrix(obs), is.array(fcst), dim(fcst)[1:2] == dim(obs))
  ## number of lead times
  nlead <- nrow(fcst)

  ## check probability inputs
  stopifnot(minprob < 0.5, minprob > 0, nprob >=1)
  ## set nodes for quantile correction
  prob <- seq(minprob, 1 - minprob, length=nprob)
  
  ## check window and jump
  stopifnot(window > 0, jump <= window)
  window <- min(window, nlead)
  
  ## check anomalies and / or debias
  if (debias & multiplicative) stop("Combination not implemented yet, do not know how to do this!")
  
  
  ## preprocess data, first de-bias if needed
  if (debias) {
    f.clim <- rowMeans(fcst, na.rm=TRUE)
    o.clim <- rowMeans(obs, na.rm=TRUE)
    if (smooth) f.clim <- sloess(f.clim, span=span)
    if (smoothobs) o.clim <- sloess(o.clim, span=span)
    if (multiplicative){
      f.clim[f.clim == 0] <- 1
      o.clim[o.clim == 0] <- 1
      f.anom <- fcst / f.clim
      o.anom <- obs / o.clim
      fout.anom <- fcst.out / f.clim
    } else {
      f.anom <- fcst - f.clim
      o.anom <- obs - o.clim
      fout.anom <- fcst.out - f.clim
    }
  } else {
    f.anom <- fcst
    o.anom <- obs
    fout.anom <- fcst.out
  }
  ## second compute signal anomalies if necessary
  if (anomalies){
    f.ens <- rowMeans(f.anom, dims=2, na.rm=T)
    fout.ens <- rowMeans(fout.anom, dims=2, na.rm=T)
    if (multiplicative){
      f.ens[f.ens == 0] <- 1
      fout.ens[fout.ens == 0] <- 1
      f.anom <- f.anom / c(f.ens)
      o.anom <- o.anom/f.ens
      fout.anom <- fout.anom / c(fout.ens)
    } else {
      f.anom <- f.anom - c(f.ens)
      o.anom <- o.anom - f.ens
      fout.anom <- fout.anom - c(fout.ens)
    }    
  }

  ## start main quantile mapping algorithm
  
  ## initialize output
  fcst.debias <- NA*fout.anom
  
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
    ## lead times for estimating quantiles
    ind <- seq(mini[i], maxi[i])
    ## lead times to be bias corrected
    ind2 <- seq(mino[i], maxo[i])

    ## compute quantiles
    oq <- quantile(o.anom[ind,], prob=prob, type=8, na.rm=T)
    fq <- quantile(f.anom[ind,,], prob=prob, type=8, na.rm=T)

    ## find boundaries in between quantiles
    fqbnds <- sort(fq[-length(fq)] + 0.5*diff(fq))
    
    ## get the probability of the output fcst given the forecast
    ## i.e. the reverse quantile function
    ## assume constant correction by discrete quantiles
    fout.qi <- findInterval(fout.anom[ind2,,], fqbnds) + 1
    if (multiplicative){
      qcorr <- oq/fq
      qcorr[fq == 0 & oq != 0] <- 1
      qcorr[fq == 0 & oq == 0] <- 0
      fcst.debias[ind2,,] <- fout.anom[ind2,,] * qcorr[fout.qi]
    } else {
      ## dry day correction
      ndry <- sum(fq == min(fq))
      if (ndry > 1){
        ## randomly assign values that fall in identical, 
        ## lowest n quantiles to quantile bins such that 
        ## these are equally populated
        ffi <- !is.na(fout.qi) & fout.qi == ndry
        ffsample <- rep(1:ndry, ceiling(sum(ffi)/ndry))
        fout.qi[ffi] <- sample(ffsample, sum(ffi), replace=FALSE)
      }
      fcst.debias[ind2,,] <- fout.anom[ind2,,] - (fq - oq)[fout.qi]    
    }      
  } ## end of loop on lead times
  
  ## postprocess output
  if (anomalies){
    if (multiplicative){
      fcst.debias <- fcst.debias * c(fout.ens)
    } else {
      fcst.debias <- fcst.debias + c(fout.ens)
    }    
  } 
  
  if (debias){
    if (multiplicative){
      fcst.debias <- fcst.debias * o.clim
    } else {
      fcst.debias <- fcst.debias + o.clim      
    }
  }
  
  if (!is.null(lower.bound)){
    fcst.debias[fcst.debias < lower.bound] <- lower.bound
  }
  
  return(fcst.debias)
}

#'@rdname qqmap
#'
qqmap_mul <- function(..., multiplicative=TRUE){
  qqmap(..., multiplicative=multiplicative)
}

#'@rdname qqmap
#'
fastqqmap <- function(fcst, window=min(nrow(fcst), 31), jump=min(window, 11), ...){
  qqmap(fcst=fcst, window=window, jump=jump, ...)
}

#'@rdname qqmap
#'
fastqqmap_mul <- function(fcst, window=min(nrow(fcst), 31), 
                          jump=min(window, 11), multiplicative=TRUE, ...){
  fastqqmap(fcst=fcst, window=window, jump=jump, multiplicative=multiplicative, ...)
}
#'@rdname qqmap
#'
qqmap_debias <- function(..., debias=TRUE){
  qqmap(..., debias=debias)
}

