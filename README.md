biascorrection
==============

[![Build Status](https://travis-ci.org/jonasbhend/biascorrection.svg?branch=master)](https://travis-ci.org/jonasbhend/biascorrection) [![codecov](https://img.shields.io/codecov/c/github/jonasbhend/biascorrection.svg)](https://codecov.io/github/jonasbhend/biascorrection)

Introduction
------------

Collection of functions to calibrate daily ensemble forecasts. This package includes a wrapper for performing bias correction in in-sample, cross-validation, moving blocks cross-validation and forward modes and also includes a variety of calibration functions (regression-based, quantile mapping). The focus is on ease-of-use rather than performance when applied to large datasets. Consequently, specific functions will have to be implemented separately to be applicable in an operational setting. No support is provided so far to apply calibration to large datasets.

The bias correction options include:

``` {.r}
library(biascorrection)
#> Loading required package: qmap
#> Loading required package: fitdistrplus
list_methods()
#>            METHODS
#> 2              ccr
#> 3      ccr_monthly
#> 4            ccrlm
#> 5                 
#> 6             comb
#> 7      conditional
#> 9      debiasApply
#> 10                
#> 11       fastqqmap
#> 12         initial
#> 13          linmod
#> 16         monthly
#> 17             mul
#> 18           qqmap
#> 20          smooth
#> 21      smooth_mul
#> 22                
#> 23    smooth_scale
#> 24       smoothccr
#> 25   smoothobs_mul
#> 26                
#> 27 smoothobs_scale
#> 28           trend
#> 29         useqmap
#>                                                  DESCRIPTION
#> 2                           Climate conserving recalibration
#> 3                  Daily CCR with monthly correction factors
#> 4                       Climate conserving recalibration (as
#> 5                                                regression)
#> 6                    Conditional bias with linear time trend
#> 7                               Bias conditional on forecast
#> 9  Apply             Apply Bias Correction to Large Ensemble
#> 10                                        Forecast Data Sets
#> 11                                          Quantile mapping
#> 12                       Bias dependent on initial condition
#> 13                                  Linear model of the bias
#> 16                       Daily calibration with monthly mean
#> 17                                 Multiplicative de-biasing
#> 18                                          Quantile mapping
#> 20                                           Mean de-biasing
#> 21                   Multiplicative de-biasing with smoothed
#> 22                                             climatologies
#> 23                     Smoothed mean de-biasing with scaling
#> 24                 Daily CCR with monthly correction factors
#> 25              Multiplicative de-biasing with smoothed obs.
#> 26                                               climatology
#> 27           Smooth mean (obs. only) de-biasing with scaling
#> 28                               Bias with linear time trend
#> 29                   Quantile mapping using the qmap package
```

Getting started
---------------

First, install the package and vignettes.

``` {.r}
devtools::install_github("jonasbhend/biascorrection", build_vignettes=TRUE)
```

Next, check out the examples provided in the vignette or the `help` and `examples` of individual functions in `biascorrection`.

``` {.r}
vignette('biascorrection')
```
