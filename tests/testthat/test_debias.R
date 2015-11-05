library(biascorrection)
context("bias correction")

fcst <- array(rnorm(1000), c(10,30,15))
fc.time <- outer(seq(1, length=nrow(fcst)), seq(1961, length=ncol(fcst)), function(x,y) as.Date(paste(y,'01', x, sep='-')))
na.vec <- rep(1, nrow(fcst))
na.vec[1] <- NA
obs <- array(rnorm(1000), c(10, 30))
fnames <- unclass(lsf.str(envir = asNamespace("biascorrection"), all=T))
mnames <- setdiff(fnames, c("debias", "list_methods", "monmean", "sloess"))
fna <- array(rnorm(1000), c(10, 30, 15))
fna[,1:10,4:15] <- NA

test_that("Missing value handling", {
  expect_error(debias(fcst, obs*na.vec))
  expect_error(debias(fcst*na.vec, obs))
  # expect_error(debias(fcst, obs, fcst.out=fcst*na.vec))
  expect_error(debias(fcst, obs, method='monthly', fc.time=fc.time*na.vec))
  })

test_that("Missing ensemble members for in-sample", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(mean(!is.na(debias(fna[,,1:3], obs, method=mn, 
                                    fc.time=fc.time, fcst.out=fna))), 
                 mean(!is.na(fna)))
  }  
})

test_that("Missing ensemble members for forward", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(mean(!is.na(debias(fna[,,1:3], obs, method=mn, 
                                    fc.time=fc.time, fcst.out=fna,
                                    forward=TRUE))), 
                 mean(!is.na(fna)))
  }  
})

test_that("Missing ensemble members for LOO crossval", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(mean(!is.na(debias(fna[,,1:3], obs, method=mn, 
                                    fc.time=fc.time, fcst.out=fna,
                                    crossval=TRUE))), 
                 mean(!is.na(fna)))
  }  
})

test_that("Missing ensemble members split sample", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(mean(!is.na(debias(fna[,,1:3], obs, method=mn, 
                                    fc.time=fc.time, fcst.out=fna,
                                    crossval=TRUE, blocklength=15))), 
                 mean(!is.na(fna)))
  }  
})

test_that('Output dimensions for in-sample', {
  for (mn in mnames){
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time)), 
                 dim(fcst))
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F])), i)
  }
})

test_that('Output dimensions for forward', {
  for (mn in mnames){
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, forward=TRUE)), 
                 dim(fcst))
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], forward=TRUE)), i)
  }
})

test_that('Output dimensions for LOO crossval', {
  for (mn in mnames){
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, crossval=TRUE)), 
                 dim(fcst))
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], crossval=TRUE)), i)
  }
})

test_that('Output dimensions for split sample', {
  for (mn in mnames){
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, crossval=TRUE, blocklength=15)), 
                 dim(fcst))
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], crossval=TRUE)), i)
  }
})
