biascorrection
==============

Introduction
------------

Collection of functions to calibrate daily ensemble forecasts. This package includes a wrapper for performing bias correction in in-sample, cross-validation, moving blocks cross-validation and forward modes and also includes a variety of calibration functions (regression-based, quantile mapping). The focus is on ease-of-use rather than performance when applied to large datasets. Consequently, specific functions will have to be implemented separately to be applicable in an operational setting. No support is provided so far to apply calibration to large datasets.

The bias correction options include:

``` {.r}
library(biascorrection)
#> Loading required package: qmap
#> Loading required package: fitdistrplus
list_methods()
#>            METHODS                                     DESCRIPTION
#> 2              ccr                Climate conserving recalibration
#> 3      ccr_monthly       Daily CCR with monthly correction factors
#> 4            ccrlm            Climate conserving recalibration (as
#> 5                                                      regression)
#> 6             comb         Conditional bias with linear time trend
#> 7      conditional                    Bias conditional on forecast
#> 9        fastqqmap                                Quantile mapping
#> 10         initial             Bias dependent on initial condition
#> 11          linmod                        Linear model of the bias
#> 14         monthly             Daily calibration with monthly mean
#> 15             mul                       Multiplicative de-biasing
#> 16           qqmap                                Quantile mapping
#> 18          smooth                                 Mean de-biasing
#> 19      smooth_mul         Multiplicative de-biasing with smoothed
#> 20                                                   climatologies
#> 21    smooth_scale           Smoothed mean de-biasing with scaling
#> 22       smoothccr       Daily CCR with monthly correction factors
#> 23   smoothobs_mul    Multiplicative de-biasing with smoothed obs.
#> 24                                                     climatology
#> 25 smoothobs_scale Smooth mean (obs. only) de-biasing with scaling
#> 26           trend                     Bias with linear time trend
#> 27         useqmap         Quantile mapping using the qmap package
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
