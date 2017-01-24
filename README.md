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
#>            METHODS                                    DESCRIPTION
#> 2              ccr               Climate Conserving Recalibration
#> 4            ccrlm           Climate Conserving Recalibration (as
#> 5                                                     Regression)
#> 6             comb        Conditional Bias With Linear Time Trend
#> 7      conditional                   Bias Conditional on Forecast
#> 11         initial            Bias dependent on initial condition
#> 12          linmod              Linear Models for Bias Correction
#> 15         monthly            Daily Calibration with Monthly Mean
#> 16             mul                      Multiplicative De-biasing
#> 17           qqmap                               Quantile Mapping
#> 19          smooth        Mean De-biasing With Smoothing of Daily
#> 20                                                    Climatology
#> 21      smooth_mul        Multiplicative De-biasing With Smoothed
#> 22                                                  Climatologies
#> 23    smooth_scale Smoothed Mean De-biasing With Variance Scaling
#> 24   smoothobs_mul        Multiplicative De-biasing With Smoothed
#> 25                                        Observation Climatology
#> 26 smoothobs_scale          Mean De-biasing With Variance Scaling
#> 27           trend                    Bias With Linear Time Trend
#> 28         useqmap      Quantile Mapping Using the 'qmap' Package
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
