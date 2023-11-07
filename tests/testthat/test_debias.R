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
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna)), 
                 !is.na(fna))
  }  
})

test_that("Missing ensemble members for forward", {
  for (mn in setdiff(mnames, c('useqmap', 'ccrlm'))){
    expect_equal(!is.na(debias(fna[,,1:3], obs, method=mn, 
                               fc.time=fc.time, fcst.out=fna,
                               strategy='forward')), 
                 !is.na(fna))
  }  
})

test_that("Missing ensemble members for LOO crossval", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(
      !is.na(
        suppressWarnings(
          debias(fna[,,1:3], obs, method=mn, 
                 fc.time=fc.time, fcst.out=fna,
                 strategy='crossval')
          
        )
      ), 
      !is.na(fna))
  }  
})

test_that("Missing ensemble members split sample", {
  for (mn in setdiff(mnames, 'useqmap')){
    expect_equal(
      !is.na(
        suppressWarnings(
          debias(fna[,,1:3], obs, method=mn, 
                 fc.time=fc.time, fcst.out=fna,
                 strategy=list(type = 'block', blocklength=15)
          )
        )
      ), 
      !is.na(fna))
  }  
})

test_that('Output dimensions for in-sample', {
  for (mn in mnames){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, method=mn, fc.time=fc.time)
    )
    expect_equal(dim(fcst_deb), dim(fcst))
    for (drop in c(TRUE, FALSE)){
      fcst_deb <- suppressWarnings(
        debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time)
      )
      expect_equal(dim(fcst_deb),
                   dim(fcst[,,1,drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        fcst_deb <- suppressWarnings(
          debias(fcst[1,,, drop=drop], obs[1,,drop=drop], method=mn, fc.time=fc.time[1,,drop=drop])
        )
        expect_equal(dim(fcst_deb), dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F])
    )
    expect_equal(ncol(fcst_deb), i)
  }
})

test_that('Output dimensions for forward', {
  for (mn in mnames){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, method=mn, fc.time=fc.time, strategy='forward')
    )
    expect_equal(dim(fcst_deb), dim(fcst))
    for (drop in c(TRUE, FALSE)){
      fcst_deb <- suppressWarnings(
        debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy='forward')        
      )
      expect_equal(dim(fcst_deb),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        fcst_deb <- suppressWarnings(
          debias(fcst[1,,, drop=drop], obs[1,,drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy='forward')
        )
        expect_equal(dim(fcst_deb),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='forward')
    )
    expect_equal(ncol(fcst_deb), i)
  }
})

test_that('Output dimensions for LOO crossval', {
  for (mn in mnames){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, method=mn, fc.time=fc.time, strategy='crossval')
    )
    expect_equal(dim(fcst_deb), 
                 dim(fcst))
    for (drop in c(TRUE, FALSE)){
      fcst_deb <- suppressWarnings(
        debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy='crossval')
      )
      expect_equal(dim(fcst_deb),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        fcst_deb <- suppressWarnings(
          debias(fcst[1,,, drop=drop], obs[1,, drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy='crossval')
        )
        expect_equal(dim(fcst_deb),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }
  for (i in 1:ncol(fcst)){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='crossval')
    )
    expect_equal(ncol(fcst_deb), i)
  }
})

test_that('Output dimensions for split sample', {
  for (mn in mnames){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, method=mn, fc.time=fc.time, strategy=list(type='block', blocklength=15))
      )
      expect_equal(dim(fcst_deb), dim(fcst))
    for (drop in c(TRUE, FALSE)){
      fcst_deb <- suppressWarnings(
        debias(fcst[,,1, drop=drop], obs, method=mn, fc.time=fc.time, strategy=list(type='block', blocklength=15))
      )
      expect_equal(dim(fcst_deb),
                   dim(fcst[,,1, drop=drop]))
      if (! mn %in% c("comb", "combRecal", "trend", "trendRecal", "conditional", "conditionalRecal")) {
        fcst_deb <- suppressWarnings(
          debias(fcst[1,,, drop=drop], obs[1,, drop=drop], method=mn, fc.time=fc.time[1,,drop=drop], strategy=list(type='block', blocklength=15))
        )
        expect_equal(dim(fcst_deb),
                     dim(fcst[1,,, drop=drop]))
      }
    }
  }    
  for (i in 1:ncol(fcst)){
    fcst_deb <- suppressWarnings(
      debias(fcst, obs, fcst.out=fcst[,1:i,,drop=F], strategy='crossval')
    )
    expect_equal(ncol(fcst_deb), i)
  }
})
