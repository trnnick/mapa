# Tests for print, summary and plot methods and tsaggr

fit <- mapaest(admissions)

test_that("summary methods return tables", {
  s <- summary(fit)
  expect_s3_class(s, "summary.mapa.fit")
  expect_equal(s$table$AL, 1:12)
  expect_true(all(c("Model","alpha") %in% colnames(s$table)))
  expect_false(any(c("nParam","AICc") %in% colnames(s$table)))
  expect_output(print(s), "ETS\\(")
  
  sc <- summary(mapacalc(admissions, fit, fh=6))
  expect_equal(dim(sc$table), c(6,5))
  expect_output(print(sc), "Combined components")
  
  sf <- summary(mapafor(admissions, fit, ifh=3, conf.lvl=0.9))
  expect_equal(rownames(sf$errors), c("t+1","t+2","t+3"))
  expect_equal(colnames(sf$errors), c("N","ME","MAE","RMSE"))
  expect_output(print(sf), "In-sample errors")
  expect_output(print(mapafor(admissions, fit, ifh=3)), "t\\+1 MSE")
})

test_that("plots run with user arguments", {
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fit, main="Test"))
  expect_silent(plot(mapaest(admissions, type="es", minimumAL=2)))
  f <- mapacalc(admissions, fit)
  expect_silent(plot(f, main="Test", xlim=c(2012,2015), col=c("black","red")))
  expect_silent(plot(f, components=TRUE))
  expect_silent(plot(f, components=c("level","season")))
  expect_error(plot(f, components="cycle"))
  fo <- mapafor(admissions, fit, ifh=3, conf.lvl=c(0.8,0.95))
  expect_silent(plot(fo, main="Test", col.PI="lightblue", insample=FALSE))
  expect_silent(plot(mapafor(admissions, fit, fh=0)))
  # User layout is respected
  par(mfrow=c(1,2))
  plot(fo)
  expect_equal(par("mfrow"), c(1,2))
  plot(f, components=TRUE)
  expect_equal(par("mfrow"), c(1,2))
})

test_that("tsaggr aggregates correctly", {
  y <- 1:13
  out <- tsaggr(y, c(2,3), fmean=FALSE)
  expect_equal(as.numeric(out$out$AL2), colSums(matrix(2:13,nrow=2)))
  expect_equal(as.numeric(out$out$AL3), colSums(matrix(2:13,nrow=3)))
  out <- tsaggr(admissions, 12, fmean=TRUE)
  expect_equal(frequency(out$out$AL12), 1)
  expect_equal(as.numeric(out$out$AL12), colMeans(matrix(tail(admissions,60),nrow=12)))
  pdf(NULL); on.exit(dev.off())
  expect_silent(tsaggr(admissions, c(2,12), outplot=TRUE, main="Test"))
})

test_that("mapafit.old gives the structure of earlier versions", {
  expect_warning(old <- mapafit.old(fit), "will be removed in a later version")
  expect_true(is.matrix(old))
  expect_s3_class(old, "mapa.fit")
  expect_equal(dim(old), c(12, 23))
  expect_equal(rownames(old), paste0("AL",1:12))
  expect_equal(colnames(old)[c(13,14,19:23)], 
               c("method","components","use","AL","original.ppy","etstype","pr.comp"))
  expect_equal(old[[1,13]], fit$models[[1]]$method)
  # The old structure is still accepted and gives the same forecasts
  expect_equal(mapacalc(admissions, old)$forecast, mapacalc(admissions, fit)$forecast)
  expect_output(print(old), "Aggregation level: 12")
  # es with skipped levels
  set.seed(1)
  x <- cbind(rnorm(74), rnorm(74))
  fitx <- mapaest(admissions, xreg=x, type="es")
  oldx <- suppressWarnings(mapafit.old(fitx))
  expect_equal(unlist(oldx[,"use"]), fitx$use, ignore_attr=TRUE)
  expect_equal(mapacalc(admissions, oldx, xreg=x)$forecast, mapacalc(admissions, fitx, xreg=x)$forecast)
})
