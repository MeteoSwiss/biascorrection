fcst <- array(1, c(1, 10, 11))
obs <- t(rnorm(10))

fcst2 <- array(rnorm(10*11), c(1,10,11))
obs2 <- t(rep(0,10))

test_that("CCR works with degenerate forecasts", {
  expect_equal(biascorrection:::ccr(fcst, obs - mean(obs), type='calibration'), fcst*0)
  expect_equal(biascorrection:::ccr(fcst, obs - mean(obs), type='prediction'), fcst*0)  
  expect_equal(biascorrection:::ccr(fcst2, obs2, type='calibration'), fcst2*0)
  expect_equal(biascorrection:::ccr(fcst2, obs2, type='prediction'), fcst2*0)  
})
