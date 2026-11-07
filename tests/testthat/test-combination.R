# Tests for the combination of components across aggregation levels

test_that("combination operators match their definitions", {
  ALs <- 1:4
  ppy <- 4
  fh <- 3
  set.seed(1)
  FCs <- array(rnorm(length(ALs)*5*fh), c(length(ALs),5,fh))
  seas <- (ppy %% ALs == 0) & (ALs < ppy)   # AL 1 and 2
  w <- function(a){(1/a)/sum(1/a)}
  wsum <- function(X,a){colSums(X*w(a))}
  
  out <- MAPA:::mapacomb(ALs,ppy,FCs,"mean")
  expect_equal(out$forecasts, colMeans(FCs[,2,]) + colMeans(FCs[,3,]) + 
                 colMeans(FCs[seas,4,]) + colMeans(FCs[,5,]))
  
  out <- MAPA:::mapacomb(ALs,ppy,FCs,"median")
  md <- function(X){apply(X,2,median)}
  expect_equal(out$forecasts, md(FCs[,2,]) + md(FCs[,3,]) + md(FCs[seas,4,]) + md(FCs[,5,]))
  
  out <- MAPA:::mapacomb(ALs,ppy,FCs,"wght")
  expect_equal(out$forecasts, wsum(FCs[,2,],ALs) + wsum(FCs[,3,],ALs) + 
                 wsum(FCs[seas,4,],ALs[seas]) + wsum(FCs[,5,],ALs))
  
  out <- MAPA:::mapacomb(ALs,ppy,FCs,"w.mean")
  expect_equal(out$forecasts, colMeans(FCs[,2,]) + colMeans(FCs[,3,]) + 
                 wsum(FCs[seas,4,],ALs[seas]) + wsum(FCs[,5,],ALs))
  
  out <- MAPA:::mapacomb(ALs,ppy,FCs,"w.median")
  expect_equal(out$forecasts, md(FCs[,2,]) + md(FCs[,3,]) + 
                 wsum(FCs[seas,4,],ALs[seas]) + wsum(FCs[,5,],ALs))
  expect_equal(colSums(out$combined), out$forecasts)
})

test_that("combination works for a single horizon and without seasonal levels", {
  FCs <- array(1:15, c(3,5,1))
  out <- MAPA:::mapacomb(5:7,1,FCs,"w.mean")
  expect_length(out$forecasts, 1)
  expect_equal(out$combined["Season",], c("t+1"=0)[[1]], ignore_attr=TRUE)
})

test_that("hybrid forecasts average the combination and the lowest level ETS", {
  fit <- mapaest(admissions)
  f1 <- mapacalc(admissions, fit, fh=6, hybrid=TRUE)
  f0 <- mapacalc(admissions, fit, fh=6, hybrid=FALSE)
  expect_equal(as.numeric(f1$forecast), 
               (f1$components[1,"ETS",] + as.numeric(f0$forecast))/2, ignore_attr=TRUE)
  expect_equal(as.numeric(f0$forecast), colSums(f0$combined), ignore_attr=TRUE)
})
