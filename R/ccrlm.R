#' Climate conserving recalibration (as regression)
#' 
#' Climate conserving recalibration with correction factors for each lead time 
#' without smoothing (application to daily series not recommended).
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param pred m vector, n x m matrix, or list of m vectors or n x m matrices
#'   with additional predictors
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param pred.out \code{n} vector, \code{n} x \code{m} matrix, or list of n
#'   vectors or n x m matrices with additional predictors
#' @param ... additional arguments for compatibility with other bias correction 
#'   methods
#'   
#' @examples
#' nveri <- 50
#' ncal <- 30
#' nn <- nveri + ncal
#' ical <- seq(1, ncal)
#' iveri <- seq(ncal + 1, nn)
#' nens <- 15
#' signal <- rnorm(nn, sd=5)
#' sdfrac <- 2
#' fcst <- array(rnorm(nn*nens, sd=sdfrac), c(1, nn, nens)) + 
#'   0.6 * signal + seq(-3,3,length=nn) 
#' obs <- array(rnorm(nn, sd=1), c(1, nn)) + signal
#' fcst.deb <- biascorrection:::ccrlm(fcst[,ical,,drop=FALSE], 
#'                                    obs[,ical,drop=FALSE],
#'                                    fcst.out=fcst)
#' fcst.deb2 <- biascorrection:::ccrlm(fcst[,ical,,drop=FALSE], 
#'                                    obs[,ical,drop=FALSE],
#'                                    fcst.out=fcst, 
#'                                    pred=seq(1,ncal), 
#'                                    pred.out=seq(1,nn))
#' f.rmse <- sqrt(apply((obs[,ical, drop=FALSE] - 
#'             apply(fcst.deb[,ical,,drop=FALSE], 2, mean))**2, 1, mean))
#' f.rmse2 <- sqrt(apply((obs[,ical,drop=FALSE] - 
#'             apply(fcst.deb2[,ical,,drop=FALSE], 2, mean))**2, 1, mean))
#' f.sd <- sqrt(mean(apply(fcst.deb[,ical,,drop=FALSE], 1:2, sd)**2))
#' f.sd2 <- sqrt(mean(apply(fcst.deb2[,ical,,drop=FALSE], 1:2, sd)**2))
#' f.sd / f.rmse ## should be one
#' f.sd2 / f.rmse2 ## should be one
#' ## spread to error with out-of-sample validation
#' f.rmse <- sqrt(apply((obs[,iveri, drop=FALSE] - 
#'             apply(fcst.deb[,iveri,,drop=FALSE], 2, mean))**2, 1, mean))
#' f.rmse2 <- sqrt(apply((obs[,iveri,drop=FALSE] - 
#'             apply(fcst.deb2[,iveri,,drop=FALSE], 2, mean))**2, 1, mean))
#' f.sd <- sqrt(mean(apply(fcst.deb[,iveri,,drop=FALSE], 1:2, sd)**2))
#' f.sd2 <- sqrt(mean(apply(fcst.deb2[,iveri,,drop=FALSE], 1:2, sd)**2))
#' f.sd / f.rmse ## should be close to one
#' f.sd2 / f.rmse2 ## can be close to one
#' 
#' @keywords util
ccrlm <- function(fcst, obs, pred=NULL, 
                  fcst.out=fcst, pred.out=pred, ...){
  warning("ccrlm is still under development, use with care")
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.spread <- sqrt(rowMeans(apply(fcst, 1:2, sd)**2))
  fcst.out.ens <- rowMeans(fcst.out, dims=2, na.rm=T)
  nfcst <- ncol(fcst)
  
  ## check input arguments and make sure it's a list with matrices
  if (!is.null(pred)){
    stopifnot(mode(pred) == mode(pred.out))
    if (is.vector(pred)){
      stopifnot(length(pred) == nfcst)
      stopifnot(length(pred.out) == ncol(fcst.out))
    } else if (is.matrix(pred)){
      stopifnot(dim(pred) == dim(fcst))
      stopifnot(dim(pred.out) == dim(fcst.out))    
    } 
    
    if (!is.list(pred)){
      pred <- list(pred=pred)
      pred.out <- list(pred=pred.out)
    }
    
    stopifnot(length(pred) == length(pred.out))
    pred <- lapply(pred, function(x) {
      if (is.vector(x)){
        x <- matrix(x, nrow(fcst), ncol(fcst), byrow=TRUE)
      }
      return(x)
    })
    pred.out <- lapply(pred.out, function(x) {
      if (is.vector(x)){
        x <- matrix(x, nrow(fcst.out), ncol(fcst.out), byrow=TRUE)
      }
      return(x)
    })
    ## make sure it is a named list
    if (is.null(names(pred))) names(pred) <- names(pred.out) <- paste(seq(pred))
    
  } ## end of if on pred
  
  
  ## simple implementation with loop
  fcst.debias <- NA*fcst.out
  for (i in 1:nrow(fcst.debias)){
    fdf <- data.frame(obs=obs[i,], fcst=fcst.ens[i,])
    odf <- data.frame(fcst=fcst.out.ens[i,])
    if (!is.null(pred)){
      for (nn in names(pred)){
        fdf[[nn]] <- pred[[nn]][i,]
        odf[[nn]] <- pred.out[[nn]][i,]
      }
    }
    flm <- lm(obs ~ ., fdf)
    pred.lm <- predict(flm, newdata=odf, interval='predict')
    tfrac <- -qt((1 - 0.95)/2, flm$df.res)
    pred.sd <- apply(pred.lm[,-1,drop=F], 1, diff) / 2 / tfrac
    repred.lm <- predict(flm, newdata=fdf, interval='predict')
    repred.sd <- apply(repred.lm[,-1], 1, diff) / 2 / tfrac
    ## fcst.debias[i,,] <- predict(flm, newdata=odf) + 
    ##   (fcst.out[i,,] - fcst.out.ens[i,]) * sd(flm$res) / fcst.spread[i] * sqrt((nfcst - 1) / nfcst)
    spread.corr <- pred.sd / fcst.spread[i]
    spread.corr[fcst.spread[i] == 0] <- 1
    fcst.debias[i,,] <- pred.lm[,1] + 
      (fcst.out[i,,] - fcst.out.ens[i,]) * spread.corr # / sqrt(mean(repred.sd**2))
  }
  
  return(fcst.debias)
}
