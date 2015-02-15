## biascorrection 0.4.7

* added lead time dependent quantile corrections in `qqmap`

## biascorrection 0.4.6

* Changed default behaviour of `debias` to allow forward calibration with new forecasts (i.e. extending past verifying observations)

* Added new CCR method based on linear regression - still under development: calibrated forecasts tend to be under-confident (in-sample)

## biascorrection 0.4.5

* Added a switch for excluding the signal (forecast ensemble mean) from quantile mapping in `qqmap`

* Debug `conditional`: prediction of bias consistently based on ensemble mean

## biascorrection 0.4.4

* Support for quantile mapping using the `qmap` package

* Added optional multiplicative quantile mapping to custom-built quantile mapping function (`qqmap`)

## biascorrection 0.4.3

* added multiplicative debiasing (e.g. for rainfall or indices measuring variability)

## biascorrection 0.4.2.1

* `debias` exclude missing value requirement for `ccr` as a fix to avoid problems in processing

## biascorrection 0.4.2

* Changed computation of spread for methods with smoothing (e.g. `smooth`, `smoothobs`). Spread is now computed as sqaure root of mean squared deviation from smoothed climatology.

* simplify missing value handling: no missing values allowed (check in `debias`)

## biascorrection 0.4.1

* added vignette on basic functionality of biascorrection package

* renamed functions: 
    * `ccr` now performs standard CCR
    * `smoothccr` and `ccr_monthly` are the corresponding routines corresponding to `monthly` and `smooth` for daily data