#' Climate conserving recalibration
#' 
#' Climate conserving recalibration with correction factors for each lead time
#' without smoothing (application to daily series not recommended).
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be
#'   applied (defaults to \code{fcst})
#' @param ... additional arguments for compatibility with other bias correction
#'   methods
#'   
#' @examples
#' nveri <- 100
#' ncal <- 10
#' nn <- nveri + ncal
#' nens <- 5
#' signal <- rnorm(nn, sd=5)
#' sdfrac <- 2
#' fcst <- array(rnorm(nn*nens, sd=sdfrac), c(1, nn, nens)) + 0.5 * signal 
#' obs <- array(rnorm(nn, sd=1), c(1, nn)) + signal
#' fcst.debias <- biascorrection:::ccrlm(fcst[,1:ncal,,drop=F], obs[,1:ncal,drop=F],
#'                                       fcst.out=fcst[, ncal + 1:nveri,,drop=F])
#' f.rmse <- sqrt(apply((obs[,ncal + 1:nveri, drop=F] - apply(fcst.debias, 2, mean))**2, 1, mean))
#' f.rmse2 <- sqrt(apply((obs[,ncal + 1:nveri,drop=F] - apply(fcst[,ncal + 1:nveri,,drop=F], 2, mean))**2, 1, mean))
#' f.sd <- sqrt(mean(apply(fcst.debias, 2, sd)**2))
#' f.sd2 <- sqrt(mean(apply(fcst[,ncal + 1:nveri,,drop=F], 2, sd)**2))
#' f.sd / f.rmse ## should be exactly 1
#' f.sd2 / f.rmse2 ## can be anything
#' mean(fcst.debias - obs[1,]) ## should be 0 (rounding errors)
#' 
#' @keywords util
ccrlm <- function(fcst, obs, fcst.out=fcst, ...){
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.spread <- sqrt(rowMeans(apply(fcst, 1:2, sd)**2))
  fcst.out.ens <- rowMeans(fcst.out, dims=2)
  
  ## simple implementation with loop
  fcst.debias <- NA*fcst.out
  for (i in 1:nrow(fcst.debias)){
    fdf <- data.frame(obs=obs[i,], fcst=fcst.ens[i,])
    flm <- lm(obs ~ fcst, fdf)
    fcst.debias[i,,] <- predict(flm, newdata=data.frame(fcst=fcst.out.ens[i,])) + 
      (fcst.out[i,,] - fcst.out.ens[i,]) * sd(flm$res) / fcst.spread[i]
  }
  
  return(fcst.debias)
  return(fi_ccrlm)
}
