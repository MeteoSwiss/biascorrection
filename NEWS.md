## biascorrection 0.5.6.9001

* added support for missing ensemble members in `fcst.out` (except for `useqmap`), forecast arrays for calibration (i.e. `fcst`) have to be complete

## biascorrection 0.5.6.9000

* changed `month` to `monmean` to avoid confusion with `month` from package `lubridate`

## biascorrection 0.5.6

* Added smooth scaling correction for both on anomalies from the climatology `smooth_scale` and anomalies from the ensemble mean `smooth_scalespread`. The latter corrects ensemble spread only with the signal being untouched.

## biascorrection 0.5.5

* Added dry-day threshold correction in `fastqqmap`. Please note, dry-day correction only works for additive corrections.

## biascorrection 0.5.4.3

* Use R's built-in weighted least squares functionality

## biascorrection 0.5.4.2

* Minor bug fix

## biascorrection 0.5.4.1

* Added option `type` also to `linmod`

## biascorrection 0.5.4

* Added option `type` to `ccr` for out-of-sample recalibration
* Fixed bug in `ccrlm` and `linmod`

## biascorrection 0.5.3.1

* Changed default for `fastqqmap` to `window = 31` and `minjump = 11`

## biascorrection 0.5.3

* Included fast quantile mapping `fastqqmap` using moving windows with larger steps

## biascorrection 0.5.2.1

* Debugged division by zero in `smooth_mul` and `smoothobs_mul`

## biascorrection 0.5.2

* correct negative climatological values for `smooth_mul` and `smoothobs_mul`

## biascorrection 0.5.1

* new differences switch for model fitting on first order differences (`linmod`)

## biascorrection 0.5

* calibration methods now use linear modelling capabilities (`linmod`)
* most calibration methods are expressed using the `linmod` function

## biascorrection 0.4.12

* new conditional approach with rescaling of residuals
* new quantile mapping with 15 day windows (doesn't work well)

## biascorrection 0.4.12

* new forward approach: first 15 years are calibrated 'backwards'

## biascorrection 0.4.11

* added CCR for daily data based on linear regression (using polynomial lead-time dependency)

## biascorrection 0.4.10

* changed treatment of first forecast in forward mode. First `nforward` forecasts are calibrated normally. Default is `nforward = 10`.

## biascorrection 0.4.9.1

* simplified notation in `conditional`, still not lead-time dependent

## biascorrection 0.4.9

* lead-time dependency of bias correction for both `qqmap` and `useqmap`
    - Changed default behaviour for `qqmap`, lead time dependency is switched of as in old version
    - Window length for lead time dependency is kept constant
    - `nn` denotes the total window length, windows with even length are supported

## biascorrection 0.4.8

* profiled `qqmap` and improved speed (less accuracy)

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