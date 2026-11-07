# Tests for estimation and forecasting

fit <- mapaest(admissions)

test_that("default horizons are consistent across functions", {
  # Default is a full season
  expect_length(mapacalc(admissions, fit)$forecast, 12)
  expect_length(mapafor(admissions, fit)$outfor, 12)
  expect_length(mapasimple(admissions)$forecast, 12)
  expect_length(mapa(admissions)$outfor, 12)
  # fh = 0 means no out-of-sample forecasts in mapafor and mapa
  f0 <- mapafor(admissions, fit, fh=0, ifh=2)
  expect_null(f0$outfor)
  expect_equal(dim(f0$infor), c(2, length(admissions)))
  expect_null(mapa(admissions, fh=0)$outfor)
  # fh = 0 with prediction intervals does not fail
  expect_null(mapafor(admissions, fit, fh=0, conf.lvl=0.9)$PI)
  # Earlier sentinels are still accepted
  expect_length(mapacalc(admissions, fit, fh=0)$forecast, 12)
  expect_length(mapafor(admissions, fit, fh=-1)$outfor, 12)
  # ppy that differs from the ts frequency gives a different default
  expect_length(mapasimple(admissions, ppy=4)$forecast, 4)
})

test_that("forecasts of ts inputs are ts and continue the series", {
  f <- mapacalc(admissions, fit, fh=5)$forecast
  expect_s3_class(f, "ts")
  expect_equal(start(f), c(2013, 6))
  expect_equal(frequency(f), 12)
  fn <- mapacalc(as.numeric(admissions), mapaest(as.numeric(admissions), ppy=12), fh=5)$forecast
  expect_false(is.ts(fn))
})

test_that("minimumAL > 1 uses the correct aggregation levels", {
  fit3 <- mapaest(admissions, minimumAL=3, maximumAL=8)
  expect_equal(fit3$AL, 3:8)
  out <- mapacalc(admissions, fit3, fh=6)
  expect_equal(dimnames(out$components)[[1]], paste0("AL",3:8))
  expect_match(capture.output(print(fit3))[2], "Aggregation level: 3")
})

test_that("skipped aggregation levels are handled", {
  # xreg reduces the sample available at high aggregation levels
  set.seed(1)
  x <- cbind(rnorm(74), rnorm(74))
  fitx <- mapaest(admissions, xreg=x, type="es")
  expect_false(all(fitx$use))
  out <- mapafor(admissions, fitx, xreg=x)
  expect_length(out$outfor, 12)
  expect_true(all(is.finite(out$outfor)))
  # The lowest aggregation level is skipped
  ys <- ts(admissions[1:26], frequency=24)
  fits <- mapaest(ys, model="ANA", type="es")
  expect_false(fits$use[1])
  expect_true(all(is.finite(mapacalc(ys, fits, fh=5)$forecast)))
})

test_that("xreg input is checked against mapafit", {
  set.seed(1)
  x <- cbind(rnorm(74), rnorm(74))
  fitx <- mapaest(admissions, xreg=x, type="es", maximumAL=6)
  expect_error(mapacalc(admissions, fitx, xreg=x[,1]), "does not match")
  expect_error(mapacalc(admissions, fitx), "provide xreg")
  expect_error(mapacalc(admissions, fitx, xreg=x[1:65,]), "must be equal")
})

test_that("principal components are the same in estimation and forecasting", {
  set.seed(1)
  n <- length(admissions)
  x <- cbind(rnorm(n+12), rnorm(n+12), rnorm(n+12))
  fitp <- mapaest(admissions, xreg=x, type="es", pr.comp=2, maximumAL=6)
  pc <- fitp$pr.comp
  # Components used in estimation, at the original frequency
  est <- MAPA:::mapaprcomp(x[1:n,], list("pr.comp"=2, "mean"=pc$mean, "sd"=pc$sd))$x.out
  for (AL in c(1,4)){
    q <- n %/% AL
    estA <- MAPA:::mapa.aggr(est, AL, q)
    forc <- MAPA:::mapaprcomp(MAPA:::mapa.aggr(x[1:n,], AL, q), pc)$x.out
    expect_equal(estA, forc, ignore_attr=TRUE)
  }
  expect_length(mapacalc(admissions, fitp, xreg=x)$forecast, 12)
})

test_that("prediction intervals are ordered", {
  out <- mapafor(admissions, fit, conf.lvl=c(0.8,0.95))
  expect_true(all(out$PI["Upper 0.95",] >= out$PI["Upper 0.80",]))
  expect_true(all(out$PI["Upper 0.80",] >= out$outfor))
  expect_true(all(out$PI["Lower 0.80",] <= out$outfor))
  expect_true(all(out$PI["Lower 0.95",] <= out$PI["Lower 0.80",]))
})

test_that("invalid options give informative errors", {
  expect_error(mapaest(admissions, paral=3), "paral")
  expect_error(mapacalc(admissions, fit, outplot=5), "outplot")
  expect_error(mapaest(as.numeric(admissions)), "ppy")
  expect_error(mapaest(admissions, paral=1), "cluster")
})

test_that("parallel processing gives the same results", {
  skip_on_cran()
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))
  fitp <- mapaest(admissions, paral=1, cl=cl)
  expect_equal(mapacalc(admissions, fitp)$forecast, mapacalc(admissions, fit)$forecast)
  expect_equal(mapafor(admissions, fit, ifh=2, paral=1, cl=cl)$infor, 
               mapafor(admissions, fit, ifh=2)$infor)
})
