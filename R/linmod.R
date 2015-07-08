#' @name linmod
#'
#' @title
#' Linear model of the bias
#' 
#' @description
#' Compute calibration for biases that follow various linear models
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble 
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param fcst.out array of forecast values to which bias correction should be 
#'   applied (defaults to \code{fcst})
#' @param fc.time forecast times as R-dates for trend aggregation
#' @param fcout.time forecast time for array to which bias correction is applied
#'   for back compatibility with leave-one-out cross-validation in
#'   \code{\link{debias}}
#' @param formula model formula to model deviations from the climatology 
#'  (see details)
#' @param recal logical, should the ensemble spread be recalibrated? (see details)
#' @param smoothobs logical, should observation climatology (and residual standard
#'   deviation) be smoothed?
#' @param smooth logical, should model climatology (and ensemble spread) be
#'   smoothed?
#' @param span the parameter which controls the degree of smoothing (see 
#'   \code{\link{loess}})
#' @param differences logical, should model be fit on first order differences (see details)?
#' @param ... additional arguments for compatibility with other calibration
#'   methods
#'   
#' @details This is the workhorse for calibration methods that can be expressed
#'   as linear models. The systematic model errors are decomposed into a
#'   seasonally varying (and thus lead-time dependent) bias and ensemble
#'   mean errors that depend on forecast anomalies from the mean forecast, on
#'   lead time, and/or on forecast date. 
#'   
#'   A variety of linear models for calibration
#'   can be specified using \code{R}'s \code{formula} notation and a few
#'   examples are given below (see also \code{example(linmod)}).
#' 
#' \describe{
#'   \item{\code{obs ~ offset(fcst) - 1}}{
#'   Simple bias correction with 
#'   seasonally varying mean bias}
#'   \item{\code{obs ~ fcst}}{
#'   Model error dependent on forecast (with 
#'   dependence being constant across lead times)}
#'   \item{\code{obs ~ fcst + fcst:poly(lead,3)}}{
#'   Model error dependent 
#'   on forecast with dependence varying across lead times}
#'   \item{\code{obs ~ offset(fcst) + year}}{
#'   Linear time trend in bias 
#'   (time trend constant across lead times)}
#'   \item{\code{obs ~ offset(fcst) + year*as.factor(format(date, "\%m"))}}{
#'   Linear time trend in bias  with trend depending on forecast month}
#' }
#' 
#' In addition to complex dependence of the systematic ensemble mean errors on 
#' forecast, lead-time, and forecast dates, \code{linmod} also allows a recalibration
#' of the ensemble spread. If \code{recal = TRUE}, the lead time dependent ensemble
#' spread is inflated (shrunk) to reflect the lead time dependent standard deviation
#' of the observation residuals. 
#' 
#' By default, the seasonally varying bias and the residual errors and model spread
#' are smoothed using a \code{loess} smoothing. Smoothing of the model ensemble mean 
#' and model spread can be switched of with \code{smooth = FALSE}, smoothing of the
#' observation climatology and residual standard deviation can be switched of with
#' \code{smoothobs = FALSE}
#'   
#' Both the residual standard deviation (\code{psd}) and the ensemble spread 
#' (\code{fsd}) are smoothed using \code{\link{loess}} as follows.
#' 
#' \code{psd.smooth = exp(loess(log(psd) ~ log(seq(along=psd)), span=span)$fit)} 
#' 
#' Alternatively, the model parameters can be estimated on time series after taking 
#' first order differences (but preserving the forecast signal) to reduce the effect
#' of auto-correlation in the residuals. 
#'   
#' @note
#' The linear models are fit using maximum likelihood assuming \code{iid} residuals
#' which will generally not apply for real forecasts (auto-correlation and
#' heteroscedasticity).
#' 
#' @examples
#' seasonal   <- sin(seq(0,4,length=215)) 
#' signal     <- outer(seasonal, rnorm(30), '+') 
#' fcst       <- array(rnorm(215*30*51), c(215,30, 15)) + 2*c(signal) 
#' obs        <- array(rnorm(215*30, mean=outer(seasonal,
#'   seq(1,3,length=30), '+')), c(215, 30)) + signal 
#' fc.time    <- outer(1:215, 1981:2010, 
#'   function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x) 
#' fdeb       <- list() 
#' fdeb$raw   <- fcst[,21:30,] ## explore calibration with varying complexity 
#' fdeb$bias  <- biascorrection:::linmod(fcst[,1:20,], obs[,1:20],
#' fcst.out=fcst[,21:30,], fc.time=fc.time[,1:20], fcout.time=fc.time[,21:30],
#' formula=obs ~ offset(fcst)) 
#' fdeb$trend <- biascorrection:::linmod(fcst[,1:20,], obs[,1:20], 
#' fcst.out=fcst[,21:30,], fc.time=fc.time[,1:20], fcout.time=fc.time[,21:30], 
#' formula=obs ~ offset(fcst) + year + year:poly(lead,3)) 
#' fdeb$cond  <- biascorrection:::linmod(fcst[,1:20,], obs[,1:20], 
#' fcst.out=fcst[,21:30,], fc.time=fc.time[,1:20], fcout.time=fc.time[,21:30], 
#' formula=obs ~ fcst + fcst:poly(lead,3)) 
#' fdeb$all   <- biascorrection:::linmod(fcst[,1:20,], obs[,1:20], 
#' fcst.out=fcst[,21:30,], fc.time=fc.time[,1:20], fcout.time=fc.time[,21:30], 
#' formula=obs ~ fcst + fcst:poly(lead,3) + year +
#' year:poly(lead,3)) 
#' fdeb$month <- biascorrection:::linmod(fcst[,1:20,], obs[,1:20], 
#' fcst.out=fcst[,21:30,], fc.time=fc.time[,1:20], fcout.time=fc.time[,21:30], 
#' formula=obs ~ fcst*as.factor(format(date, '%m')))
#' 
#' fdeb.rmse  <- lapply(fdeb, function(x) sqrt(apply((rowMeans(x,dims=2) -
#' obs[,21:30])**2, 1, mean))) 
#' boxplot(fdeb.rmse, ylab='RMSE of the ensemble
#' mean')
#' 
#' @keywords util
#' @export
linmod <- function(fcst, obs, fcst.out=fcst, 
                   fc.time=NULL,
                   fcout.time=NULL,
                   formula=obs ~ fcst, 
                   recal=FALSE, 
                   smooth=TRUE, smoothobs=TRUE,
                   span=min(1, 91/nrow(fcst)), 
                   differences=FALSE, ...){
  
  ## internal function
  fdate <- function(x,y){
    as.Date(paste0(y,'-01-01')) + x - 1
  }
  
  if (is.null(fc.time)) fc.time <- outer(1:nrow(fcst), 1980 + 1:ncol(fcst), fdate)
  if (is.null(fcout.time)) fcout.time <- outer(1:nrow(fcst.out), 1980 + 1:ncol(fcst.out), fdate)
  
  stopifnot(length(fcout.time) == length(fcst.out[,,1]))
  stopifnot(length(fc.time) == length(obs))
  stopifnot(dim(fcst)[1:2] == dim(obs))
  fcst.ens <- rowMeans(fcst, dims=2)
  fcst.ens[is.na(obs)] <- NA
  fcst.mn <- rowMeans(fcst.ens, dims=1, na.rm=T)
  fcst.clim <- if (smooth) sloess(fcst.mn, span=span) else fcst.mn
  obs.mn <- rowMeans(obs, dims=1, na.rm=T)
  obs.clim <- if (smoothobs) sloess(obs.mn, span=span) else obs.mn
  
  fcst.out.ens <- rowMeans(fcst.out, dims=2)
  
  in.df <- data.frame(fcst=c(fcst.ens - fcst.clim),
                      obs=c(obs - obs.clim), 
                      lead=seq(1, nrow(obs)),
                      year=as.numeric(format(fc.time[1:nrow(obs), 1:ncol(obs)], '%Y')),
                      date=fc.time[1:nrow(obs), 1:ncol(obs)])
  out.df <- data.frame(fcst=c(fcst.out.ens - fcst.clim),
                       lead=seq(1,nrow(fcst.out)),
                       year=as.numeric(format(fcout.time, '%Y')),
                       date=fcout.time)
  
  f.lm <- lm(formula, in.df)

  if (differences){
    ## estimate first-order autoregression of residuals
    fres <- array(f.lm$res, dim(obs))
    pp <- mean(apply(fres, 2, function(x) pacf(x, plot=F)$acf[1]))
    
    ## set up new input data frame
    in.df2 <- data.frame(fcst=c(fcst.ens[-1,] - pp*fcst.ens[-nrow(fcst.ens),] + 
                                  rep(colMeans(fcst.ens), each=nrow(fcst.ens)-1)),
                         obs=c(obs[-1,] - pp*obs[-nrow(obs),] + 
                                 rep(colMeans(obs), each=nrow(obs) - 1)),
                         lead=seq(1, nrow(obs) - 1),
                         year=as.numeric(format(fc.time[-1,], '%Y')),
                         date=fc.time[-1,])
    f.lm <- lm(formula, in.df2)
  }
  
  ## compute lead-time dependent inflation for recalibration
  if (recal){
    fres <- obs - predict(f.lm, newdata=in.df)
    psd <- apply(fres, 1, sd)
    fsd <- apply(fcst.out - c(fcst.out.ens), 1, sd)
    if (smooth) fsd <- exp(loess(log(fsd) ~ log(seq(fsd)), span=span)$fit)
    if (smoothobs) psd <- exp(loess(log(psd) ~ log(seq(psd)), span=span)$fit)
    inflate <- psd / fsd 
      
  } else {
    inflate <- 1
  }
  
  fcst.debias <- obs.clim + 
    predict(f.lm, newdata=out.df) + 
    (fcst.out - c(fcst.out.ens)) * inflate
  
  return(fcst.debias)
}
