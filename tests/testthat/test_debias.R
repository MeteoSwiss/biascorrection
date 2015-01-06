library(biascorrection)
context("bias correction")

fcst <- array(rnorm(1000), c(10,30,51))
fc.time <- outer(seq(1, length=nrow(fcst)), seq(1961, length=ncol(fcst)), function(x,y) as.Date(paste(y,'01', x, sep='-')))
na.vec <- rep(1, nrow(fcst))
na.vec[1] <- NA
obs <- array(rnorm(1000), c(10, 30))
mnames <- setdiff(ls(pos='package:biascorrection'), c('debias', 'month', 'sloess', 'cor'))

test_that("Missing value handling", {
  expect_error(debias(fcst, obs*na.vec))
  expect_error(debias(fcst*na.vec, obs))
  expect_error(debias(fcst, obs, fcst.out=fcst*na.vec))
  expect_error(debias(fcst, obs, method='monthly', fc.time=fc.time*na.vec))
  })

test_that('Output dimensions are correct', {
  for (mn in mnames){
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time)), dim(fcst))
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F])), i)
  }
})


