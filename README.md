biascorrection
==============

[![Build Status](https://travis-ci.org/jonasbhend/biascorrection.svg?branch=master)](https://travis-ci.org/jonasbhend/biascorrection) [![codecov](https://img.shields.io/codecov/c/github/jonasbhend/biascorrection.svg)](https://codecov.io/github/jonasbhend/biascorrection)

Introduction
------------

Collection of functions to calibrate daily ensemble forecasts. This package includes a wrapper for performing bias correction in in-sample, cross-validation, moving blocks cross-validation and forward modes and also includes a variety of calibration functions (regression-based, quantile mapping). The focus is on ease-of-use rather than performance when applied to large datasets. Consequently, specific functions will have to be implemented separately to be applicable in an operational setting.

Marginal support is provided to apply calibration functions to large datasets. Function `debiasApply` automates calibration with collections of forecasts and calibrating observations.

The bias correction options include:

``` r
library(biascorrection)
#> Loading required package: qmap
#> Loading required package: fitdistrplus
#> Loading required package: easyVerification
#> Loading required package: SpecsVerification
list_methods()
#>  METHODS          DESCRIPTION                                   
#>  ccr              Climate Conserving Recalibration              
#>  ccrlm            Climate Conserving Recalibration (as          
#>                   Regression)                                   
#>  comb             Conditional Bias With Linear Time Trend       
#>  conditional      Bias Conditional on Forecast                  
#>  initial          Bias dependent on initial condition           
#>  linmod           Linear Models for Bias Correction             
#>  monthly          Daily Calibration with Monthly Mean           
#>  moving           Moving Window Mean De-biasing                 
#>  mul              Multiplicative De-biasing                     
#>  qqmap            Quantile Mapping                              
#>  smooth           Mean De-biasing With Smoothing of Daily       
#>                   Climatology                                   
#>  smooth_mul       Multiplicative De-biasing With Smoothed       
#>                   Climatologies                                 
#>  smooth_scale     Smoothed Mean De-biasing With Variance Scaling
#>  smoothobs_mul    Multiplicative De-biasing With Smoothed       
#>                   Observation Climatology                       
#>  smoothobs_scale  Mean De-biasing With Variance Scaling         
#>  trend            Bias With Linear Time Trend                   
#>  useqmap          Quantile Mapping Using the 'qmap' Package
```

Getting started
---------------

First, install the package and vignettes.

``` r
devtools::install_github("jonasbhend/biascorrection", build_vignettes=TRUE)
```

Next, check out the examples provided in the vignette or the `help` and `examples` of individual functions in `biascorrection`.

``` r
vignette('biascorrection')
```
