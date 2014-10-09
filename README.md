biascorrection
==============

Collection of functions to de-bias daily ensemble forecasts. This package includes a variety of functions with a focus on ease-of-use rather than performance when applied to large datasets. Accordingly, specific functions will have to be implemented separately to be applicable in an operational setting.

The bias correction options include:
* mean de-biasing (unbias)
* mean de-biasing with loess smoothing (smooth, smoothobs)
* monthly-mean de-biasing (monthly)
* climate conserving recalibration based on smoothed daily climatology (ccr)
* climate conserving recalibration based on monthly climatology (ccr_monthly)
* parametric quantile mapping based on smoothed climatology and variance (qqmap)
* semi-parametric quantile mapping (qqmap_semi)
* mean de-biasing with linear time trend (trend)
