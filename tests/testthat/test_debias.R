library(biascorrection)
context("bias correction")

fcst <- array(rnorm(1000), c(10,30,15))
obs <- array(rnorm(1000), c(10, 30))
fc.time <- outer(seq(1, length=nrow(fcst)), seq(1961, length=ncol(fcst)), function(x,y) as.Date(paste(y,'01', x, sep='-')))
na.vec <- rep(1, nrow(fcst))
na.vec[1] <- NA
na.vec2 <- rep(1, length(obs))
na.vec2[rep(c(1:4, rep(NA, 26)), each=nrow(fcst))] <- NA ## throws error if less than 5 forecasts
fnames <- unclass(lsf.str(envir = asNamespace("biascorrection"), all=T))
mnames <- setdiff(fnames, c("debias", "list_methods", "monmean", "sloess", "iqqmap", "print.biascorrection"))
fna <- array(rnorm(1000), c(10, 30, 15))
fna[,1:10,4:15] <- NA

test_that("Missing ensemble members for in-sample", {
  for (mn in setdiff(mnames, c('useqmap', 'ccrlm'))){
    print(mn)
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna)), 
                 !is.na(fna))
  }  
})

test_that("Missing ensemble members for forward", {
  for (mn in setdiff(mnames, c('useqmap', 'ccrlm'))){
    print(mn)
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna,
                               strategy='forward')), 
                 !is.na(fna))
  }  
})

test_that("Missing ensemble members for LOO crossval", {
  for (mn in setdiff(mnames, 'useqmap')){
    print(mn)
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna,
                               strategy='crossval')), 
                 !is.na(fna))
  }  
})

test_that("Missing ensemble members split sample", {
  for (mn in setdiff(mnames, 'useqmap')){
    print(mn)
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna,
                               strategy=list(type = 'block', blocklength=15))), 
                 !is.na(fna))
  }  
})

test_that('Output dimensions for in-sample', {
  for (mn in mnames){
    print(mn)
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time)), 
                 dim(fcst))
    for (drop in c(TRUE, FALSE)){
      expect_equal(dim(debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time)),
                   dim(fcst[,,1,drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        expect_equal(dim(debias(fcst[1,,, drop=drop], obs[1,,drop=drop], method=mn, fc.time=fc.time[1,,drop=drop])),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F])), i)
  }
})

test_that('Output dimensions for forward', {
  for (mn in mnames){
    print(mn)
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, strategy='forward')), 
                 dim(fcst))
    for (drop in c(TRUE, FALSE)){
      expect_equal(dim(debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy='forward')),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        expect_equal(dim(debias(fcst[1,,, drop=drop], obs[1,,drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy='forward')),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='forward')), i)
  }
})

test_that('Output dimensions for LOO crossval', {
  for (mn in mnames){
    print(mn)
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, strategy='crossval')), 
                 dim(fcst))
    for (drop in c(TRUE, FALSE)){
      expect_equal(dim(debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy='crossval')),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        expect_equal(dim(debias(fcst[1,,, drop=drop], obs[1,, drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy='crossval')),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='crossval')), i)
  }
})

test_that('Output dimensions for split sample', {
  for (mn in mnames){
    print(mn)
    expect_equal(dim(debias(fcst, obs, method=mn, fc.time=fc.time, strategy=list(type='block', blocklength=15))), 
                 dim(fcst))
    for (drop in c(TRUE, FALSE)){
      expect_equal(dim(debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy=list(type='block', blocklength=15))),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        expect_equal(dim(debias(fcst[1,,, drop=drop], obs[1,, drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy=list(type='block', blocklength=15))),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }    
  for (i in 1:ncol(fcst)){
    expect_equal(ncol(debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='crossval')), i)
  }
})
