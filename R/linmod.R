#' @name linmod
#'   
#' @title Linear model of the bias
#'   
#' @description Compute calibration for biases that follow various linear models
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
#' @param formula model formula to model deviations from the climatology (see 
#'   details)
#' @param recal logical, should the ensemble spread be recalibrated? (see 
#'   details)
#' @param smoothobs logical, should observation climatology (and residual 
#'   standard deviation) be smoothed?
#' @param smooth logical, should model climatology (and ensemble spread) be 
#'   smoothed?
#' @param span the parameter which controls the degree of smoothing (see 
#'   \code{\link{loess}})
#' @param bleach logical, should variance of the residuals be controlled for 
#'   (see details)?
#' @param differences logical, should model be fit on first order differences 
#'   (see details)?
#' @param type one of \code{prediction} or \code{calibration}, where signal 
#'   calibration uncertainties are propagated for \code{type = "prediction"} 
#'   (see details).
#' @param ... additional arguments for compatibility with other calibration 
#'   methods
#'   
#' @details This is the workhorse for calibration methods that can be expressed 
#'   as linear models. The systematic model errors are decomposed into a 
#'   seasonally varying (and thus lead-time dependent) bias and ensemble mean 
#'   errors that depend on forecast anomalies from the mean forecast, on lead 
#'   time, and/or on forecast date.
#'   
#'   A variety of linear models for calibration can be specified using 
#'   \code{R}'s \code{formula} notation and a few examples are given below (see 
#'   also \code{example(linmod)}).
#'   
#'   \describe{ \item{\code{obs ~ offset(fcst) - 1}}{ Simple bias correction 
#'   with seasonally varying mean bias} \item{\code{obs ~ fcst}}{ Model error 
#'   dependent on forecast (with dependence being constant across lead times)} 
#'   \item{\code{obs ~ fcst + fcst:poly(lead,3)}}{ Model error dependent on 
#'   forecast with dependence varying across lead times} \item{\code{obs ~ 
#'   offset(fcst) + year}}{ Linear time trend in bias (time trend constant 
#'   across lead times)} \item{\code{obs ~ offset(fcst) + 
#'   year*as.factor(format(date, "\%m"))}}{ Linear time trend in bias  with 
#'   trend depending on forecast month} }
#'   
#'   In addition to complex dependence of the systematic ensemble mean errors on
#'   forecast, lead-time, and forecast dates, \code{linmod} also allows a 
#'   recalibration of the ensemble spread. If \code{recal = TRUE}, the lead time
#'   dependent ensemble spread is inflated (shrunk) to reflect the lead time 
#'   dependent standard deviation of the observation residuals. In addition, if 
#'   \code{type = 'prediction'}, the ensemble spread matches the predictive 
#'   standard deviation for the out-of-sample forecast. That is, for a simple 
#'   linear model with \eqn{y = ax + b + \varepsilon}{y = ax + b + e}, the 
#'   predictive standard error is \eqn{\sigma \sqrt{1 + 1/n + {x_0}^2 / 
#'   \sum{x_j^2}}}{ s * sqrt(1 + 1/n + x0^2 / sum(xj^2))}, where \eqn{n} is the 
#'   number of forecasts in the calibration set, \eqn{x_0}{x0} is the ensemble 
#'   mean of the foracst that is to be bias-corrected, and \eqn{x_j}{xj} are the
#'   forecast ensemble means in the calibration set.
#'   
#'   The underlying assumption of \code{iid} residuals in the linear model is 
#'   generally not fullfilled, thus weighted least squares can be used instead 
#'   with \code{bleach = TRUE}. The weighting is chosen proportional to the
#'   inverse of the lead-time dependent variance of the residuals.
#'   
#'   By default, the seasonally varying bias and the residual errors and model 
#'   spread are smoothed using a \code{loess} smoothing. Smoothing of the model 
#'   ensemble mean and model spread can be switched of with \code{smooth = 
#'   FALSE}, smoothing of the observation climatology and residual standard 
#'   deviation can be switched of with \code{smoothobs = FALSE}
#'   
#'   Both the residual standard deviation (\code{psd}) and the ensemble spread 
#'   (\code{fsd}) are smoothed using \code{\link{loess}} as follows.
#'   
#'   \code{psd.smooth = exp(loess(log(psd) ~ log(seq(along=psd)), 
#'   span=span)$fit)}
#'   
#'   Alternatively, the model parameters can be estimated on time series after 
#'   taking first-order differences (but preserving the forecast signal) to 
#'   reduce the effect of auto-correlation in the residuals with 
#'   \code{differences = TRUE}.
#'   
#' @note The linear models are fit using maximum likelihood assuming \code{iid} 
#'   residuals which will generally not apply for real forecasts 
#'   (auto-correlation and heteroscedasticity).
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
                   bleach=TRUE, 
                   differences=FALSE, 
                   type=c("calibration", "prediction"), 
                   ...){
  
  type <- match.arg(type)
  
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
    pp <- tanh(mean(atanh(apply(fres, 2, function(x) pacf(x, plot=F)$acf[1]))))
    
    ## use lead-time deviation from seasonal mean signal
    ffres <- fcst.ens - fcst.clim - rep(colMeans(fcst.ens - fcst.clim), each=nrow(fcst.ens)) 
    oores <- obs - obs.clim - rep(colMeans(obs - obs.clim), each=nrow(obs))
    
    in.df2 <- data.frame(fcst=c(ffres[-1,] - pp*ffres[-nrow(ffres),] + 
                                  rep(colMeans(fcst.ens - fcst.clim), each=nrow(ffres) - 1)),
                         obs=c(oores[-1,] - pp*oores[-nrow(oores),] + 
                                 rep(colMeans(obs - obs.clim), each=nrow(oores) - 1)),
                         lead=seq(1, nrow(obs) - 1),
                         year=as.numeric(format(fc.time[-1,], '%Y')),
                         date=fc.time[-1,])
    f.lm <- lm(formula, in.df2)
    
    if (bleach){
      sd.res <- apply(array(in.df2$obs - predict(f.lm, newdata=in.df2), dim(obs[-1,])), 1, sd)
      if (smooth) sd.res <- exp(loess(log(sd.res) ~ log(seq(sd.res)))$fit)
      in.df2$ww <- 1 / sd.res**2
      sd.res <- sd.res[c(seq(sd.res), length(sd.res))]
      f.lm <- lm(formula, in.df2, weights=ww)      
    } else {
      sd.res <- 1
    }    
  } else {
    ## compute pre-whitening and rescale forecast and obs anomalies
    if (bleach){
      sd.res <- apply(array(f.lm$res, dim(obs)), 1, sd)
      if (smooth) sd.res <- exp(loess(log(sd.res) ~ log(seq(sd.res)))$fit)
      in.df$ww <- 1 / sd.res**2
      f.lm <- lm(formula, in.df, weights=ww)
    } else {
      sd.res <- rep(1, nrow(obs))
    }
    
  }
  
  ## compute lead-time dependent inflation for recalibration
  if (recal){
    fsd <- apply(fcst - c(fcst.ens), 1, sd)
    if (smooth) fsd <- exp(loess(log(fsd) ~ log(seq(fsd)), span=span)$fit)
    if (type == 'prediction'){
      ## prediction interval is tfrac*sd_pred
      plm <- predict(f.lm, newdata=out.df, interval='prediction', level=pnorm(1), weights=1 / sd.res**2)
      tfrac <- -qt((1 - pnorm(1))/2, f.lm$df.residual)
      psd <- array((plm[,'upr'] - plm[,'fit'])/tfrac, dim(fcst.out.ens))
      if (differences){
        ## additional correction to take into account that
        ## sd(f.lm$res) != sd(in.df$obs - predict(f.lm, in.df))
        sd.corr1 <- apply(matrix(in.df$obs - predict(f.lm, in.df, weights=1/sd.res**2), nrow(obs)), 1, sd)
        sd.corr2 <- apply(matrix(f.lm$res, nrow(obs) - 1), 1, sd)
        sd.corr2 <- sd.corr2[c(seq(sd.corr2), length(sd.corr2))]
        sd.corr <- sd.corr1 / sd.corr2
        psd <- psd*sd.corr
      }
      if (smoothobs) psd <- c(apply(psd, 2, function(y) exp(loess(log(y) ~ log(seq(y)), span=span)$fit)))
    } else {
      fres <- array(in.df$obs - predict(f.lm, newdata=in.df, weights=1 / sd.res**2), dim(obs))
      psd <- apply(fres, 1, sd)
      if (smoothobs) psd <- exp(loess(log(psd) ~ log(seq(psd)), span=span)$fit)
    }
    inflate <- c(psd / fsd) 
  } else {
    inflate <- 1
  }
  
  fcst.debias <- obs.clim + 
    predict(f.lm, newdata=out.df, weights=1 / sd.res**2) + 
    (fcst.out - c(fcst.out.ens)) * inflate
  
  return(fcst.debias)
}
