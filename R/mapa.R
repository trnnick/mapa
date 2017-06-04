# Nikolaos Kourentzes & Fotios Petropoulos (2014)
# Now split into multiple files (mapa.R, mapa.estim.R, mapa.forecast.R)

#-------------------------------------------------
mapa <- function(y, ppy=NULL, fh=ppy, ifh=1, minimumAL=1, maximumAL=ppy, 
	comb=c("mean","median","w.mean","w.median","wght"), paral=c(0,1,2), display=c(0,1), outplot=c(1,0), 
	hybrid=c(TRUE,FALSE), model="ZZZ", type=c("ets","es"), conf.lvl=NULL, 
	xreg=NULL, pr.comp=0, ...){
# Wrapper to estimate and produce MAPA in- and out-of-sample forecasts
# Uses mapaest and mapafor
#  
# Inputs:
#   y           = In sample observations of a time series (vector)
#                 If y == "paper" then it prints paper reference
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If y is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   fh          = Forecast horizon. Default = ppy
#   ifh         = In-sample forecast horizon. Default = 1
#   minimumAL   = Lowest aggregation level to use. Default = 1
#   maximumAL   = Highest aggregation level to use. Default = ppy, maximumAL>1
#   comb        = Combination operator. This can be: "mean"; "median"; "wght" - where each 
#                 aggregation level is weighted inversly to aggregation; "w.mean" - level 
#                 and trend components are averaged, but seasonal and xreg follow the wght 
#                 combination; "w.median" - as w.mean, but with median. It is suggested that 
#                 for data with high sampling frequency to use one of the "w.mean" and "w.median".
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
#   model       = Allow only that type of ETS at each aggregation level. This follows similar
#                 coding to the ets function. The first letter refers to the error 
#                 type ("A", "M" or "Z"); the second letter refers to the trend 
#                 type ("N","A","Ad","M","Md", "Z", "X" or "Y"); and the third letter 
#                 refers to the season type ("N","A","M", "Z", "X" or "Y"). The letters mean: 
#                 "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, 
#                 "X"=automatically select between none and additive and "Y"=automatically
#                 select between none and multiplicative. "X" and "Y" supported only by 
#                 type=="es". If used with type=="ets" a warning will be given and they will 
#                 default to "Z".  A "d" for trend implies damped. 
#                 By default model="ZZZ". If due to sample limitation ETS cannot be calculated 
#                 at an aggregation level for the selected model, then no estimation is done 
#                 for that specific level. For aggregation levels that seasonality becomes 
#                 1 then a non-seasonal model is estimated.
#   type        = What type of exponential smoothing implementation to use. "es" - use from 
#                 the smooth package; "ets" use from the forecast package.
#   conf.lvl    = Vector of confidence level for prediction intervals. Values must be (0,1). 
#                 If conf.lvl == NULL then no intervals are calculated. For example to get 
#                 the intervals for 80% and 95% use conf.lvl=c(0.8,0.95).
#   xreg        = Vector or matrix of exogenous variables to be included in the MAPA. 
#                 If matrix then rows are observations and columns are variables. 
#                 Must be at least as long as in-sample. Additional observations are unused.
#                 Note that including xreg will force type="es". 
#   pr.comp     = MAPAx can use principal component analysis to preprocess xreg. When comp is -1 
#                 then the number of retained components is chosen automatically. When comp=0 
#                 then no pre-processing is performed and the original xreg is used. Any other 
#                 value represents the number of principal components retained. 
#   ...         = Pass additional arguments to es and ets.
#
# Output:
#   out$infor   = In-sample forecasts
#   out$outfor  = Out-of-sample forecasts
#   out$PI      = Prediction intervals for given confidence levels
#   out$MSE     = In-sample MSE error
#   out$MAE     = In-sample MAE error
  
  # Defaults
  comb <- match.arg(comb,c("mean","median","w.mean","w.median","wght"))
  paral <- paral[1]
  display <- display[1]
  outplot <- outplot[1]
  hybrid <- hybrid[1]
  type <- type[1]
  
  # Paper info
  if (!is.numeric(y)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
                     " 30 (2), 291-302.",sep=""))
    return(invisible())
  }
  
  # Estimate MAPA
  mapafit <- mapaest(y=y, ppy=ppy, minimumAL=minimumAL, 
                     maximumAL=maximumAL, paral=paral, display=display, 
                     outplot=outplot, model=model, type=type, xreg=xreg, pr.comp=pr.comp, ...)
  
  # If not horizon was given forecasts a full season
  if (is.null(fh)){
    idx.ppy <- which(colnames(mapafit)=="original.ppy")
    fh <- as.numeric(mapafit[1,idx.ppy])
  }
  
  # Produce in- and out-of-sample forecasts
  out <- mapafor(y=y, mapafit=mapafit, fh=fh, ifh=ifh, comb=comb, 
                 outplot=outplot, hybrid=hybrid, conf.lvl=conf.lvl, xreg=xreg)
  
  return(out)
  
}

#-------------------------------------------------
mapacomb <- function(minimumAL,maximumAL,ppy,FCs,comb){
# This function combines the translated ets states
  
  # perm_levels is not needed for forecasting. This is already checked in the estimation.
  # perm_levels <- array(0, maximumAL) # permitted levels due to ETS implementation (observations>=4)
  perm_seas <- array(0, maximumAL)   # permitted seasonalities
  for (AL in minimumAL:maximumAL){
    # if (observations %/% AL >=4){
    #   perm_levels[AL] <- 1
    # }
    if ((ppy %% AL == 0) & (AL<ppy)){
      perm_seas[AL] <- 1
    }
  }
  # perm_levels <- perm_levels[minimumAL:maximumAL]
  perm_levels <- rep(1,(maximumAL-minimumAL+1))
  perm_seas <- perm_seas[minimumAL:maximumAL]
  
  if (dim(FCs)[3] != 1){ # Forecast multiple steps ahead
    level <- FCs[perm_levels==1, 2, ]
    trend <- FCs[perm_levels==1, 3, ]
    season <- FCs[(perm_levels==1 & perm_seas==1), 4, ]
    xreg <- FCs[perm_levels==1, 5, ]
    # Check that all are arrays
    if (!is.array(level)){
      level <- array(level,c(1,length(level)))
    }
    if (!is.array(trend)){
      trend <- array(trend,c(1,length(trend)))
    }
    if (!is.array(season)){
      season <- array(season,c(1,length(season)))
    }
    if (!is.array(xreg)){
      xreg <- array(xreg,c(1,length(xreg)))
    }
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(colMeans(level),colMeans(trend),
                                 colMeans(season),colMeans(xreg)),na.rm=TRUE) # MAPA(mean) forecasts
    } else if (comb=="median"){
      # forecasts <- colSums(rbind(colMedians(level),colMedians(trend),
      #                            colMedians(season)),na.rm=TRUE) # MAPA(median) forecasts
      forecasts <- colSums(rbind(apply(level,2,"median"),apply(trend,2,"median"),
                                 apply(season,2,"median"),apply(xreg,2,"median")),na.rm=TRUE) # MAPA(median) forecasts
    } else if (comb=="wght"){
      # Weighted sum
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      fh <- dim(level)[2]
      wghts.level <- matrix(rep(wghts.level,fh),ncol=fh)
      wghts.season <- matrix(rep(wghts.season,fh),ncol=fh)
      forecasts <- colSums(rbind(colSums(level * wghts.level),colSums(trend * wghts.level),
                                 colSums(season * wghts.season),colSums(xreg * wghts.level)),na.rm=TRUE)
    } else if (comb=="w.mean"){
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      fh <- dim(level)[2]
      wghts.season <- matrix(rep(wghts.season,fh),ncol=fh)
      forecasts <- colSums(rbind(colMeans(level),colMeans(trend),
                                 colSums(season * wghts.season),colSums(xreg * wghts.level)),na.rm=TRUE)
    } else if (comb=="w.median"){
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      fh <- dim(level)[2]
      wghts.season <- matrix(rep(wghts.season,fh),ncol=fh)
      forecasts <- colSums(rbind(apply(level,2,"median"),apply(trend,2,"median"),
                                 colSums(season * wghts.season),colSums(xreg * wghts.level)),na.rm=TRUE)
    }
  } else {
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(mean(FCs[perm_levels==1, 2, ]),mean(FCs[perm_levels==1, 3, ]),
		    mean(FCs[(perm_levels==1 & perm_seas==1), 4, ]),mean(FCs[perm_levels==1, 5, ])), na.rm=TRUE) # MAPA(mean) forecasts
    } else if (comb=="median"){
      forecasts <- colSums(rbind(median(FCs[perm_levels==1, 2, ]),median(FCs[perm_levels==1, 3, ]),
		    median(FCs[(perm_levels==1 & perm_seas==1), 4, ]),median(FCs[perm_levels==1, 5, ])), na.rm=TRUE) # MAPA(median) forecasts
    } else if (comb=="wght"){
      # Weighted sum
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      forecasts <- sum(c(sum(FCs[perm_levels==1, 2, ]*wghts.level), sum(FCs[perm_levels==1, 3, ]*wghts.level),
                         sum(FCs[(perm_levels==1 & perm_seas==1), 4, ]*wghts.season),
                         sum(FCs[perm_levels==1, 5, ]*wghts.level)),na.rm=TRUE)
    } else if (comb=="w.mean"){
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      forecasts <- sum(c(mean(FCs[perm_levels==1, 2, ]),mean(FCs[perm_levels==1, 3, ]),
                         sum(FCs[(perm_levels==1 & perm_seas==1), 4, ]*wghts.season),
                         sum(FCs[perm_levels==1, 5, ]*wghts.level)),na.rm=TRUE)
    } else if (comb=="w.median"){
      wghts <- 1/(minimumAL:maximumAL)
      wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
      wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
      forecasts <- sum(c(median(FCs[perm_levels==1, 2, ]),median(FCs[perm_levels==1, 3, ]),
                         sum(FCs[(perm_levels==1 & perm_seas==1), 4, ]*wghts.season),
                         sum(FCs[perm_levels==1, 5, ]*wghts.level)),na.rm=TRUE)
    }
  }
  
  # Return output
  return(list(forecasts=forecasts,perm_levels=perm_levels,perm_seas=perm_seas))
}

#-------------------------------------------------
mapaplot <- function(outplot,FCs,minimumAL,maximumAL,perm_levels,perm_seas,
                     observations,y,forecasts,fh,comb){
# Produce MAPA forecast & components plot 
# outplot == 0, no plot; == 1 series plot; == 2 component plot  

  if (outplot > 0){
    FClevel <- FCs[perm_levels==1, 2, ]
    FCtrend <- FCs[perm_levels==1, 3, ]
    if (sum(perm_levels==1 & perm_seas==1)!=0){
      FCseason <- FCs[(perm_levels==1 & perm_seas==1), 4, ]
    } else {
      FCseason <- NULL
    }
    FCxreg <- FCs[perm_levels==1, 5, ]
    if (all(FCxreg == 0)){
      FCxreg <- NULL
    }

    # Check that all are arrays
    if (!is.array(FClevel)){
      FClevel <- array(FClevel,c(length(FClevel)/fh,fh))
    }
    if (!is.array(FCtrend)){
      FCtrend <- array(FCtrend,c(length(FCtrend)/fh,fh))
    }
    if (!is.array(FCseason) && !is.null(FCseason)){
      FCseason <- array(FCseason,c(length(FCseason)/fh,fh))
    }
    if (!is.array(FCxreg) && !is.null(FCxreg)){
      FCxreg <- array(FCxreg,c(length(FCxreg)/fh,fh))
    }
    
    clrs <- colorRampPalette(brewer.pal(11,"Spectral")[c(1:5,8:11)])(length(perm_levels))
    if (outplot == 2){
      if (!is.null(FCseason) & !is.null(FCxreg)){
        layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = TRUE))
      } else if (!is.null(FCseason) | !is.null(FCxreg)){
        layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
      } else {
        layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      }
    } else {
      layout(matrix(1, 1, 1, byrow = TRUE))
    }
    # Find min max
    ymax <- max(c(max(forecasts),max(y)))
    ymin <- min(c(min(forecasts),min(y)))
    yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
    # Plot prediction
    cmp <- brewer.pal(3,"Set1")
    plot(1:observations,y, type="l", col=cmp[2], xlab="", ylab="", 
		  main="Forecast", xlim=c(1, observations+fh), ylim=yminmax)
    lines((observations):(observations+fh),c(y[observations],forecasts), col=cmp[1])
    if (outplot == 2){
      # Plot Level
      ymin <- min(FClevel)
      ymax <- max(FClevel)
      yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
      plot(FClevel[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Level", ylim=yminmax)
      if (sum(perm_levels)>1){
        for (i in 2:sum(perm_levels)){
          lines(FClevel[i, ], type="l", col=clrs[i])
        }
      }
      if (comb=="mean"){
        FClevel.plot <- colMeans(FClevel)
      } else if (comb=="median"){
        FClevel.plot <- apply(FClevel,2,"median")
      } else {
        wghts <- 1/(minimumAL:maximumAL)
        wghts.level <- wghts[perm_levels==1]/sum(wghts[perm_levels==1])
        wghts.level <- matrix(rep(wghts.level,fh),ncol=fh)
        FClevel.plot <- colSums(FClevel * wghts.level)
      }
      lines(FClevel.plot, type="l", col="black", lwd=2)
      # Plot trend
      ymin <- min(FCtrend)
      ymax <- max(FCtrend)
      yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
      plot(FCtrend[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Trend", 
		    ylim=yminmax)
      if (sum(perm_levels)>1){
        for (i in 2:sum(perm_levels)){
          lines(FCtrend[i, ], type="l", col=clrs[i])
        }
      }
      if (comb=="mean"){
        FCtrend.plot <- colMeans(FCtrend)
      } else if (comb=="median"){
        FCtrend.plot <- apply(FCtrend,2,"median")
      } else {
        FCtrend.plot <- colSums(FCtrend * wghts.level)
      }
      lines(FCtrend.plot, type="l", col="black", lwd=2)
      # Plot season
      if (!is.null(FCseason)){
        ymin <- min(FCseason)
        ymax <- max(FCseason)
        yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
        plot(FCseason[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Season", ylim=yminmax)
        if (dim(FCseason)[1]>1){
          for (i in 2:sum(perm_levels==1 & perm_seas==1)){
            lines(FCseason[i, ], type="l", col=clrs[i])
          }
        }
        if (comb=="mean"){
          FCseason.plot <- colMeans(FCseason)
        } else if (comb=="median"){
          FCseason.plot <- apply(FCseason,2,"median")
        } else {
          wghts.season <- wghts[perm_levels==1 & perm_seas==1]/sum(wghts[perm_levels==1 & perm_seas==1])
          wghts.season <- matrix(rep(wghts.season,fh),ncol=fh)
          FCseason.plot <- colSums(FCseason * wghts.season)
        }
        lines(FCseason.plot, type="l", col="black", lwd=2)
      }
      # Plot xreg
      if (!is.null(FCxreg)){
        ymin <- min(FCxreg)
        ymax <- max(FCxreg)
        yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
        plot(FCxreg[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Xreg", ylim=yminmax)
        if (dim(FCxreg)[1]>1){
          for (i in 2:sum(perm_levels==1)){
            lines(FCxreg[i, ], type="l", col=clrs[i])
          }
        }
        if (comb=="mean"){
          FCxreg.plot <- colMeans(FCxreg)
        } else if (comb=="median"){
          FCxreg.plot <- apply(FCxreg,2,"median")
        } else {
          FCxreg.plot <- colSums(FCxreg * wghts.level)
        }
        lines(FCxreg.plot, type="l", col="black", lwd=2)
      }
    }
  }
}

#-------------------------------------------------
mapasimple <- function(y, ppy=NULL, fh=ppy, minimumAL=1, maximumAL=ppy, comb=c("mean","median","w.mean","w.median","wght"), 
                       paral=c(0,1,2), display=c(0,1), outplot=c(1,0), hybrid=c(TRUE,FALSE), 
                       model="ZZZ", type=c("ets","es"), xreg=NULL, pr.comp=0, ...){
  # MAPA estimation and forecast
  #  
  # Inputs:
  #   y           = In sample observations of a time series (vector)
  #                 If y == "paper" then it prints paper reference
  #   ppy         = Periods in a season of the time series at the sampled frequency.
  #                 If y is a ts object then this is taken from its frequency,
  #                 unless overriden. 
  #   fh          = Forecast horizon. Default = ppy
  #   minimumAL   = Lowest aggregation level to use. Default = 1, maximumAL>1
  #   maximumAL   = Highest aggregation level to use. Default = ppy
  #   comb        = Combination operator. This can be: "mean"; "median"; "wght" - where each 
  #                 aggregation level is weighted inversly to aggregation; "w.mean" - level 
  #                 and trend components are averaged, but seasonal and xreg follow the wght 
  #                 combination; "w.median" - as w.mean, but with median. It is suggested that 
  #                 for data with high sampling frequency to use one of the "w.mean" and "w.median".
  #   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
  #                 2 = yes and initialise cluster. Default is 0.
  #   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
  #   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
  #                 2 = time series, forecasts and components. For the components the rainbow 
  #                 colouring scheme is used. Red is aggregation level 1, followed by yellow, 
  #                 green, cyan, blue and magenta for the higher aggregation levels. Default is 1. 
  #   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE  
  #                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
  #   model       = Allow only that type of ETS at each aggregation level. This follows similar
  #                 coding to the ets function. The first letter refers to the error 
  #                 type ("A", "M" or "Z"); the second letter refers to the trend 
  #                 type ("N","A","Ad","M","Md", "Z", "X" or "Y"); and the third letter 
  #                 refers to the season type ("N","A","M", "Z", "X" or "Y"). The letters mean: 
  #                 "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, 
  #                 "X"=automatically select between none and additive and "Y"=automatically
  #                 select between none and multiplicative. "X" and "Y" supported only by 
  #                 type=="es". If used with type=="ets" a warning will be given and they will 
  #                 default to "Z".  A "d" for trend implies damped. 
  #                 By default model="ZZZ". If due to sample limitation ETS cannot be calculated 
  #                 at an aggregation level for the selected model, then no estimation is done 
  #                 for that specific level. For aggregation levels that seasonality becomes 
  #                 1 then a non-seasonal model is estimated.
  #   type        = What type of exponential smoothing implementation to use. "es" - use from 
  #                 the smooth package; "ets" use from the forecast package
  #   xreg        = Vector or matrix of exogenous variables to be included in the MAPA. 
  #                 If matrix then rows are observations and columns are variables. 
  #                 Must be at least as long as in-sample. Additional observations are unused.
  #                 Note that including xreg will force type="es". 
  #   pr.comp     = MAPAx can use principal component analysis to preprocess xreg. When comp is -1 
  #                 then the number of retained components is chosen automatically. When comp=0 
  #                 then no pre-processing is performed and the original xreg is used. Any other 
  #                 value represents the number of principal components retained. 
  #   ...         = Pass additional arguments to es and ets.
  #
  # Output:
  #   forecasts   = Vector with forecasts
  #   components  = array with MAPA components
  
  # Paper info
  if (!is.numeric(y)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
                     " 30 (2), 291-302.",sep=""))
    return(invisible())
  }  
  
  # Defaults
  comb <- match.arg(comb,c("mean","median","w.mean","w.median","wght"))
  paral <- paral[1]
  display <- display[1]
  outplot <- outplot[1]
  hybrid <- hybrid[1]
  type <- type[1]
  
  # Estimate MAPA
  mapafit <- mapaest(y=y, ppy=ppy, minimumAL=minimumAL, 
                     maximumAL=maximumAL, paral=paral, display=display, 
                     outplot=outplot, model=model, type=type, xreg=xreg, pr.comp=pr.comp, ...)
  
  # If not horizon was given forecasts a full season
  if (is.null(fh)){
    idx.ppy <- which(colnames(mapafit)=="original.ppy")
    fh <- as.numeric(mapafit[1,idx.ppy])
  }
  
  # Produce in- and out-of-sample forecasts
  out <- mapacalc(y=y, mapafit=mapafit, fh=fh, comb=comb, 
                  outplot=outplot, hybrid=hybrid, xreg=xreg)
  
  return(out)
  
}

#-------------------------------------------------
mapaprcomp <- function(x,pr.comp){
  # Function to preprocess xreg with principal component
  #
  # Inputs
  #   x          = xreg
  #   pr.comp    = list containing list("pr.comp", "mean", "sd")
  
  p <- dim(x)[2]
  
  # Get bits from list
  comp <- pr.comp$pr.comp
  x.mn <- pr.comp$mean
  x.sd <- pr.comp$sd
  
  if (comp != 0){
    # Reduction is perfomed
    
    # Normalise xreg
    n <- dim(x)[1]
    x <- (x - matrix(rep(x.mn,n),nrow=n,byrow=TRUE))/matrix(rep(x.sd,n),nrow=n,byrow=TRUE)
    
    # Check if requested number of components is possible
    if (comp > p){
      stop("The requested number of principal components exceeds the number of xreg variables.")
    }
    
    # Calculate principal components
    x.pr <- prcomp(x,center=FALSE,scale.=FALSE)
    
    # If a nmber of components is provided then simply choose that and proceed
    if (comp < 0){
      # Select number of components automatically
      # For this I am using Jolliffe's rule 
      # Cangelosi, Richard, and Alain Goriely. "Component retention in principal component analysis with application to cDNA microarray data." Biology direct 2.1 (2007): 1.
      ev <- x.pr$sdev^2
      ev <- (1-1/nrow(x))*ev
      comp <- sum(ev >= 0.7*mean(ev))
    }
    x.out <- x.pr$x[,1:comp,drop=FALSE]
    
  } else {
    # No pre-processing
    x.out <- x
  }
  
  return(list("x.out"=x.out,"pr.comp"=list("pr.comp"=comp,"mean"=x.mn,"sd"=x.sd)))
  
}