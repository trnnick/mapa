# These functions handle forecasting using MAPA
# Nikolaos Kourentzes 2016

statetranslate <- function(fit,AL,fh,q,ppyA,ets.type,xreg=NULL){
  # This function prepares ets states for MAPA combination
  # It extrapolates from last states the forecasts and translates to additive
  
  FCs_temp <- array(0, c(5, fh))
  fhA <- ceiling(fh/AL)   # forecast horizon for the aggregated level
  
  # Get model components
  components <- mapa.components(fit,ets.type)
  
  # Get last states, so that the same calculations can be used for both types
  if (ets.type == "ets"){
    # Handle ets from forecast package
  
    param <- fit$par
    pnames <- names(param)
    if (sum(pnames=="phi")>0) {
      phi <- param[pnames=="phi"]
    } else {
      phi <- NULL
    }
    
    level <- fit$states[q+1, 1]
    if (components[2]=="N"){ 
      b = 0 # indicates that there is no trend
    } else {
      b <- 1 # indicates that there is trend
      trend <- fit$states[q+1, 2]
    }
    if (components[3]!="N"){
      season <- rev(fit$states[q+1,(2+b):(ppyA+1+b)])
    }
    
  } else if (ets.type == "es"){
    # Handle es from smooth package
    init.obs <- max(1,length(fit$initial$seasonal))
    phi <- fit$phi
    
    # components <- substring(components,5,nchar(components))
    # components <- substring(components,1,nchar(components)-1)
    
    # Estimates for the Level Component - Careful es updates future level!!!
    level <- fit$states[q+init.obs, 1]
    if (components[2]!="N"){
      trend <- fit$states[q+init.obs, 2]
    }
    # sc <- dim(fit$states)[2]
    sc <- which(colnames(fit$states) == "seasonal")
    if (components[3]!="N"){
      season <- fit$states[(q+init.obs-ppyA+1):(q+init.obs),sc]
    }
    
  }
  
  # Estimates for the Level Component
  FCs_temp[2, ] <- as.numeric(rep(rep(level, fhA), each=AL)[1:fh])
  
  # Estimates for the Trend Component
  if (components[2]=="A"){ # additive trend
    FCs_temp[3, ] <- as.numeric(rep(trend * (1:fhA), each=AL))[1:fh]
  } else if (components[2]=="Ad"){ # additive damped trend
    # We divide with phi because of an internal calculation for 
    # the damped trend in the ETS package
    FCs_temp[3, ] <- as.numeric(rep(cumsum((trend/phi)* phi^(1:fhA)), each=AL))[1:fh]
  } else if (components[2]=="M"){ # multiplicative trend
    FCs_temp[3, ] <- as.numeric(rep((trend^(1:fhA)-1), each=AL)[1:fh] * FCs_temp[2,])
  } else if (components[2]=="Md"){ # multiplicative damped trend
    # We divide with phi because of an internal calculation for the damped trend in the ETS package
    FCs_temp[3, ] <- as.numeric(rep((((trend ^ (1/phi)) ^ cumsum(phi^(1:fhA)))-1), 
                                    each=AL)[1:fh] * FCs_temp[2, ])
  } else { # no trend
    # No trend
    FCs_temp[3, ] <- 0
  }
  
  # Estimates for the Seasonal Component 
  if (components[3]=="N"){ # no seasonality
    FCs_temp[4, ] <- 0
  } else if (components[3]=="A"){ # additive seasonality
    FCs_temp[4, ] <- as.numeric(rep(rep(season, fhA), each=AL))[1:fh]
  } else { # multiplicative seasonality
    FCs_temp[4, ] <- as.numeric((rep(rep(season, fhA),
                                     each=AL)[1:fh] - 1)) * (FCs_temp[2, ] + FCs_temp[3, ])
  }
  
  # Estimate for the xreg 
  # No longer paer of es - added externally
  # if (dim(fit$data)[2]>1){
  #   xreg <- fit$data[,-1,drop=FALSE]  
  # } else {
  #   xreg <- NULL
  # }
  
  if (ets.type == "es" && !is.null(xreg)){
    # Check sample once the es() data is correctly carried forward
    if (!is.null(dim(xreg))){
      FCs_temp[5, ] <- as.numeric(rep(colSums(t(xreg * matrix(rep(fit$initial$xreg,fhA),nrow=fhA,byrow=TRUE))) , each=AL))[1:fh]
      # FCs_temp[5, ] <- as.numeric(rep(colSums(t(xreg[(q+1):(q+fhA),] * matrix(rep(fit$initial$xreg,fhA),nrow=fhA,byrow=TRUE))) , each=AL))[1:fh]
    } else {
      FCs_temp[5, ] <- as.numeric(rep(xreg * fit$initial$xreg, each=AL))[1:fh]
      # FCs_temp[5, ] <- as.numeric(rep(xreg[(q+1):(q+fhA)] * fit$initial$xreg, each=AL))[1:fh]
    }
  }
  
  # Recreate ETS forecasts
  FCs_temp[1, ] <- colSums(FCs_temp[2:5,,drop=FALSE])     
  
  # Return output
  return(FCs_temp)
  
}

#-------------------------------------------------
#' MAPA Internal Calculation of Forecasts
#'
#' @description
#' The following function produces a simple MAPA forecast from a given origin.
#' This is meant to be an internal function. Use mapafor instead.
#'
#' @param y In sample observations of a time series (vector).
#' @param mapafit Fitted MAPA model (from mapaest).
#' @param fh Forecast horizon. Default = ppy. For compatibility with earlier versions, fh = 0 also gives the default.
#' @param comb Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
#' @param outplot Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, forecasts, and components. For the components the spectral colouring scheme is used. Red is the lowest aggregation level and blue the highest. Default is 0.
#' @param hybrid Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
#' @param xreg Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample plus fh. Additional observations are unused.
#'
#' @return
#' \item{forecast}{Vector with forecasts. If y is a ts object, then this is also a ts object.}
#' \item{components}{Array with MAPA components.}
#' \item{combined}{Combined components across aggregation levels.}
#' The output is of class \code{mapa.calc}, which also stores the time series and the settings used. It has \code{print}, \code{summary}, and \code{plot} methods, see \code{\link{plot.mapa.calc}}.
#'
#' @details
#' If \code{mapafit} was estimated with xreg, then xreg must also be provided here, with the same variables and at least length(y) + fh observations, so that the values for the forecast horizon are available. See \code{\link{mapaest}} for details on the use of xreg.
#'
#' \code{mapafit} can also have the structure of earlier versions (<= 2.0.6), for example a fitted MAPA saved with an earlier version, or produced by \code{\link{mapafit.old}}. Note that \code{mapafit.old} is provided only for compatibility and will be removed in a later version of the package.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}; Fotios Petropoulos.
#'
#' @seealso \code{\link{mapafor}}, \code{\link{mapa}}.
#'
#' @keywords ~mapafor
#' @examples
#' mapafit <- mapaest(admissions)
#'
#' # Forecast a full season ahead (default)
#' out <- mapacalc(admissions, mapafit)
#' print(out)
#'
#' # Forecast 18 steps ahead, and plot the forecasts and the components
#' out <- mapacalc(admissions, mapafit, fh=18, outplot=2)
#' summary(out)
#' plot(out, components=c("level","season"), main="Admissions")
#'
#' # Different combination operators, with and without the hybrid forecasts
#' out.mean <- mapacalc(admissions, mapafit, comb="mean")
#' out.wght <- mapacalc(admissions, mapafit, comb="wght", hybrid=FALSE)
#'
#' # The same works for MAPA estimated with es from the smooth package
#' mapafit.es <- mapaest(admissions, type="es")
#' out.es <- mapacalc(admissions, mapafit.es, fh=18)
#'
#' # MAPAx: xreg is needed both in estimation and forecasting. For forecasting 
#' # it must contain the in-sample observations and the fh future values.
#' set.seed(1)
#' fh <- 12
#' x <- rnorm(length(admissions)+fh)
#' mapafit.x <- mapaest(admissions, type="es", xreg=x)
#' out.x <- mapacalc(admissions, mapafit.x, fh=fh, xreg=x)
#' plot(out.x, components=TRUE)
#' @export
mapacalc <- function(y, mapafit, fh=NULL, comb=c("w.mean", "w.median", "mean", "median", "wght"),
                     outplot=c(0,1,2), hybrid=c(TRUE,FALSE), xreg=NULL){
  comb <- match.arg(comb,c("w.mean","w.median","mean","median","wght"))
  outplot <- mapa.choice(outplot,c(0,1,2),"outplot")
  hybrid <- mapa.choice(hybrid,c(TRUE,FALSE),"hybrid")
  
  # Find locations in mapafit
  # Convert mapafit from earlier versions
  mapafit <- mapafit.convert(mapafit)
  
  # Get settings from mapafit
  # Only the estimated aggregation levels are used
  idx.use <- which(mapafit$use)
  ALs <- mapafit$AL[idx.use]
  ppy <- mapafit$ppy
  ets.type <- mapafit$type
  
  # Set default foreast horizon
  if (is.null(fh) || fh == 0){
    fh <- ppy
  }
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Convert xreg to array, check its size and mapafit compatibility
  if (!is.null(xreg)){
    # Check if model is es
    if (ets.type != "es"){
      stop('Only mapafit estimated with type=="es" accepts xreg inputs.')
    }
    # Force xreg to be a matrix, each column a variable
    xreg <- mapa.xreg(xreg)
    p <- dim(xreg)[2]
    # Check that it contains at least n observatoins and trim
    xn <- dim(xreg)[1]
    if (xn < (fh + observations)){
      stop("Number of observations of xreg must be equal of exceed the length of y + fh.")
    } 
    # Check that mapafit expected xreg matches xreg input
    # Either pr.comp = 0 and check number of x's or pr.comp > 0 and x's should be >= pr.comp
    # The number of original xreg variables is now stored in mapafit, so this is a single check
    if (p != mapafit$xreg.n){
      stop("Number of xreg input variables does not match mapafit specification.")
    }
  } else if (mapafit$xreg.n > 0){
    stop("mapafit was estimated with xreg. Please provide xreg.")
  }
  
  # The forecasted components are saved here
  FCs <- array(0, c(length(ALs), 5, fh),dimnames=list(paste("AL",ALs,sep=""),
                                                      c("ETS","Level","Trend","Season","Xreg"),paste("t+",1:fh,sep=""))) 
  
  # MAPA forecast
  for (ALi in 1:length(ALs)){
    
    AL <- ALs[ALi]
    AL.fit <- mapafit$models[[idx.use[ALi]]]
    fhA <- ceiling(fh/AL)   # This is used for es and xreg preprocessing
    alinfo <- mapa.alinfo(AL,observations,ppy)
    q <- alinfo$q
    ppyA <- alinfo$ppyA
    
    # Aggregation
    yA <- mapa.aggr(y,AL,q)
    ats <- ts(yA, frequency = ppyA) 
    
    # Aggregate xreg 
    # This is aligned with the end point of the aggregation of y
    if (!is.null(xreg)){
      # Observations discarded from the beginning, so that the aggregation buckets 
      # end at the last in-sample observation. The future buckets are shifted by 
      # the same r, so that each one contains AL observations.
      r <- observations - q*AL
      temp <- xreg[(r+1):min(c(r+(q+fhA)*AL,xn)),,drop=FALSE]
      m <- (q+fhA)*AL - dim(temp)[1]
      temp <- rbind(temp,matrix(NA,m,p))
      xregA <- apply(temp,2,function(xk){colMeans(matrix(xk,nrow=AL),na.rm=TRUE)})
      xregA <- matrix(xregA,nrow=q+fhA)

      # Preprocess if needed
      pr.out <- mapaprcomp(xregA,mapafit$pr.comp)
      xregA <- pr.out$x.out
      
    } else {
      xregA <- NULL
    }
    
    # Predict
    if (ets.type == "ets"){
      
      ats.fit <- ets(ats, AL.fit, use.initial.values=TRUE)
      xregF <- NULL

    } else if (ets.type=="es"){
      
      ats.fit <- suppressWarnings(es(ats, model=AL.fit, h=fhA, silent=TRUE,xreg=xregA))
      # Smooth no longer store the future xreg, pass it as a separate object
      xregF <- tail(xregA,fhA)
      
    }
    
    # Transalte ets states for MAPA
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,ets.type,xreg=xregF)
    
    # Return MAPA components
    FCs[ALi, , ] <- FCs_temp
    
  }
  
  # # Get fit errors
  # if (ets.type == "ets"){
  #   fit.error <- unlist(mapafit[,colnames(mapafit)=="sigma2"])
  # } else {
  #   fit.error <- unlist(mapafit[,colnames(mapafit)=="s2"])
  # }
  
  # MAPA combination
  combres <- mapacomb(ALs,ppy,FCs,comb)
  forecasts <- combres$forecasts
  
  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }  
  
  # Construct output
  out <- structure(list("forecast"=mapa.tsout(forecasts,y),"components"=FCs,
                        "combined"=combres$combined,"y"=y,"fh"=fh,"comb"=comb,
                        "hybrid"=hybrid,"AL"=ALs,"ppy"=ppy,"type"=ets.type,
                        "perm_seas"=combres$perm_seas),class="mapa.calc")
  
  # Plot output
  # outplot == 0, no plot; == 1 series plot; == 2 component plot
  if (outplot > 0){
    plot(out,components=(outplot==2))
  }
  
  return(out)
  
}

#-------------------------------------------------
#' @rdname summary.mapa.fit
#' @export
print.mapa.calc <- function(x,...){
  
  cat(paste0("MAPA forecasts (",x$type,", ",x$comb,ifelse(x$hybrid,", hybrid",""),")\n"))
  fc <- x$forecast
  if (!isa(fc,"ts")){
    names(fc) <- paste("t+",1:length(fc),sep="")
  }
  print(fc)
  
  return(invisible(x))
  
}

#-------------------------------------------------
#' @rdname summary.mapa.fit
#' @export
summary.mapa.calc <- function(object,...){
  tbl <- t(rbind(object$combined,"Forecast"=as.numeric(object$forecast)))
  rownames(tbl) <- paste("t+",1:object$fh,sep="")
  
  out <- structure(list("table"=tbl,"comb"=object$comb,"hybrid"=object$hybrid,
                        "AL"=object$AL,"type"=object$type,"ppy"=object$ppy,"fh"=object$fh),
                   class="summary.mapa.calc")
  
  return(out)
  
}

#-------------------------------------------------
#' @export
print.summary.mapa.calc <- function(x,digits=max(3,getOption("digits")-3),...){
  
  cat(paste0("MAPA forecasts using ",x$type,"\n"))
  cat(paste0("Aggregation levels: ",paste(x$AL,collapse=", "),"\n"))
  cat(paste0("Combination: ",x$comb,"\tHybrid: ",x$hybrid,"\tHorizon: ",x$fh,"\n\n"))
  cat("Combined components and forecasts:\n")
  tbl <- x$table
  if (all(tbl[,"Xreg"]==0)){
    tbl <- tbl[,colnames(tbl)!="Xreg",drop=FALSE]
  }
  print(round(tbl,digits))
  
  return(invisible(x))
  
}

#-------------------------------------------------
#' Plot MAPA forecasts and components
#'
#' @description The following function plots the forecasts produced by mapacalc or mapasimple, and optionally the forecasted components at each aggregation level. The black line in each component panel is the combination of the component across aggregation levels, as used for the forecasts. The plot is the same as the one that can be produced optionally by mapacalc.
#'
#' @param x MAPA forecasts (from mapacalc or mapasimple).
#' @param components Plot the components. FALSE = no; TRUE = all available components. Alternatively a vector with any of "level", "trend", "season", and "xreg". Components that are not available are not plotted.
#' @param col Vector of 2 colours for the time series and the forecasts. If NULL the default colours are used.
#' @param col.AL Vector of colours for the components at each aggregation level. If NULL the spectral colouring scheme is used, with red the lowest and blue the highest aggregation level.
#' @param ... Additional arguments passed to \code{plot} for the forecast panel, for example \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim}, and \code{lwd}.
#'
#' @note If y is a ts object, its time index is used for the horizontal axis. When only the forecasts are plotted the current graphical layout is respected, so the plot can be placed in a user defined grid. When components are plotted the layout is restored on exit.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @seealso \code{\link{mapacalc}}, \code{\link{mapasimple}}.
#'
#' @keywords ~plot
#' @examples
#' mapafit <- mapaest(admissions)
#' out <- mapacalc(admissions, mapafit, fh=18)
#' plot(out)
#' plot(out, components=TRUE)
#' plot(out, components="season", main="Admissions", col=c("black","red"))
#'
#' # Custom colours for the aggregation levels
#' plot(out, components=c("level","trend"), col.AL=rainbow(12))
#'
#' # Place the forecast plot in a user defined layout
#' par(mfrow=c(1,2))
#' plot(out, xlim=c(2012,2015), main="Zoomed")
#' plot(mapafit)
#' par(mfrow=c(1,1))
#' @export
plot.mapa.calc <- function(x,components=c(FALSE,TRUE),col=NULL,col.AL=NULL,...){
  y <- x$y
  observations <- length(y)
  forecasts <- as.numeric(x$forecast)
  fh <- x$fh
  FCs <- x$components
  nAL <- length(x$AL)
  
  # Find which components to plot
  comps.all <- c("Level","Trend","Season","Xreg")
  if (is.logical(components)){
    if (components[1]){
      comps <- comps.all
    } else {
      comps <- NULL
    }
  } else {
    comps <- comps.all[match(tolower(components),tolower(comps.all))]
    if (any(is.na(comps))){
      stop('components must be TRUE, FALSE, or any of "level", "trend", "season", "xreg".')
    }
  }
  # Drop components that are not available 
  if (!any(x$perm_seas)){
    comps <- setdiff(comps,"Season")
  }
  if (all(FCs[,"Xreg",] == 0)){
    comps <- setdiff(comps,"Xreg")
  }
  
  # Colours
  if (is.null(col)){
    col <- brewer.pal(3,"Set1")[c(2,1)]
  }
  if (is.null(col.AL)){
    col.AL <- colorRampPalette(brewer.pal(11,"Spectral")[c(1:5,8:11)])(nAL)
  }
  col.AL <- rep(col.AL,length.out=nAL)
  
  # Layout: forecast on top, components below
  nc <- length(comps)
  if (nc > 0){
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    layout(matrix(c(rep(1,nc),2:(nc+1)), 2, nc, byrow = TRUE))
  }
  
  # Plot prediction
  tm <- mapa.time(y,fh)
  args <- list(x=tm$tin, y=as.numeric(y), type="l", col=col[1], xlab="", ylab="", 
               main="Forecast", xlim=range(c(tm$tin,tm$tout)), 
               ylim=mapa.ylim(y,forecasts))
  args <- modifyList(args,list(...))
  do.call(plot,args)
  lwd <- if (is.null(args$lwd)) 1 else args$lwd
  lines(c(tm$tin[observations],tm$tout),c(y[observations],forecasts), col=col[2], lwd=lwd)
  
  # Plot components
  # Plot Level
  # Plot trend
  # Plot season
  # Plot xreg
  # These are done in a single loop, as the code is the same for all components
  for (cmp in comps){
    if (cmp == "Season"){
      idx <- which(x$perm_seas)
    } else {
      idx <- 1:nAL
    }
    FCk <- matrix(FCs[idx,cmp,],nrow=length(idx),ncol=fh)
    plot(FCk[1, ], type="l", col=col.AL[idx[1]], xlab="", ylab="", main=cmp, 
         ylim=mapa.ylim(FCk,x$combined[cmp,]))
    if (length(idx)>1){
      for (i in 2:length(idx)){
        lines(FCk[i, ], type="l", col=col.AL[idx[i]])
      }
    }
    lines(x$combined[cmp,], type="l", col="black", lwd=2)
  }
  
  return(invisible())
  
}

#-------------------------------------------------
#' MAPA Forecast
#'
#' @description
#' The following function produces in- and out-of-sample MAPA forecasts, for multiple steps ahead.
#'  This is the recommended function to use in forecasting with MAPA.
#'
#' @param y In sample observations of a time series (vector).
#' @param mapafit Fitted MAPA model (from mapaest).
#' @param fh Forecast horizon. Default = ppy. If fh = 0 then only in-sample forecasts are produced. For compatibility with earlier versions, fh = -1 also gives the default.
#' @param ifh In-sample forecast horizon. Default = 1.
#' @param comb Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
#' @param outplot Provide output plot. 0 = no; 1 = yes. Default is 0.
#' @param hybrid Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
#' @param conf.lvl Vector of confidence level for prediction intervals. Values must be (0,1). If conf.lvl == NULL then no intervals are calculated. For example to get the intervals for 80\% and 95\% use conf.lvl=c(0.8,0.95).
#' @param xreg Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample plus fh. Additional observations are unused.
#' @param paral Use parallel processing for the in-sample forecasts. 0 = no; 1 = yes (requires initialised cluster); 2 = yes and initialise cluster. Default is 0.
#' @param cl Cluster to use when paral=1, created with \code{parallel::makeCluster}. If NULL, a cluster named \code{cl} in the global environment is used, for compatibility with earlier versions.
#'
#' @return
#' \item{infor}{In-sample forecasts.}
#' \item{outfor}{Out-of-sample forecasts.}
#' \item{PI}{Prediction intervals for given confidence levels.}
#' \item{MSE}{In-sample MSE error.}
#' \item{MAE}{In-sample MAE error.}
#' If y is a ts object, then outfor is also a ts object. The output is of class \code{mapa.for}, which also stores the time series and the settings used. It has \code{print}, \code{summary}, and \code{plot} methods, see \code{\link{plot.mapa.for}}. The summary reports in-sample errors for each horizon.
#'
#' @details
#' If \code{mapafit} was estimated with xreg, then xreg must also be provided here, with the same variables and at least length(y) + fh observations, so that the values for the forecast horizon are available. See \code{\link{mapaest}} for details on the use of xreg.
#'
#' \code{mapafit} can also have the structure of earlier versions (<= 2.0.6), for example a fitted MAPA saved with an earlier version, or produced by \code{\link{mapafit.old}}. Note that \code{mapafit.old} is provided only for compatibility and will be removed in a later version of the package.
#'
#' @note The calculation of the prediction intervals is based on the empirical multiple step ahead MSE. To speed up calculations set \code{conf.lvl=NULL}. If very long forecast horizons are requested then once no more t+h MSE can be calculated the following approximation is used: \eqn{\sqrt{MSE_{t+1}}*\sqrt{h}}{sqrt(MSE_[t+1])*sqrt(h)} for the error.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}; Fotios Petropoulos.
#'
#' @seealso \code{\link{mapa}}, \code{\link{mapaest}}, \code{\link{mapacalc}}.
#'
#' @keywords ~mapa ~mapacalc ~mapaest
#' @examples
#' mapafit <- mapaest(admissions)
#'
#' # Out-of-sample forecasts a full season ahead and in-sample 1-step ahead forecasts
#' out <- mapafor(admissions, mapafit)
#'
#' # Forecast 18 steps ahead, with in-sample forecasts up to 3 steps ahead and
#' # 80% and 95% prediction intervals
#' out <- mapafor(admissions, mapafit, fh=18, ifh=3, conf.lvl=c(0.8,0.95))
#' print(out)
#' summary(out)
#' plot(out, main="Admissions", col.PI="lightblue")
#'
#' # Only in-sample forecasts, for example to evaluate the fit up to 12 steps ahead
#' out <- mapafor(admissions, mapafit, fh=0, ifh=12)
#' summary(out)
#'
#' \donttest{
#' # With es from the smooth package and the w.median combination
#' mapafit.es <- mapaest(admissions, type="es")
#' out.es <- mapafor(admissions, mapafit.es, comb="w.median")
#'
#' # MAPAx: xreg must contain the in-sample observations and the fh future values
#' set.seed(1)
#' fh <- 12
#' x <- rnorm(length(admissions)+fh)
#' mapafit.x <- mapaest(admissions, type="es", xreg=x)
#' out.x <- mapafor(admissions, mapafit.x, fh=fh, xreg=x, conf.lvl=0.9)
#' plot(out.x)
#'
#' # Parallel processing for the in-sample forecasts, with an existing cluster
#' library(parallel)
#' cl <- makeCluster(2)
#' out <- mapafor(admissions, mapafit, ifh=3, paral=1, cl=cl)
#' stopCluster(cl)
#' }
#' @export
mapafor <- function(y, mapafit, fh=NULL, ifh=1, comb=c("w.mean","w.median","mean","median","wght"), 
                    outplot=c(0,1), hybrid=c(TRUE,FALSE), conf.lvl=NULL, xreg=NULL, 
                    paral=c(0,1,2), cl=NULL) {
  comb <- match.arg(comb,c("w.mean","w.median","mean","median","wght"))
  outplot <- mapa.choice(outplot,c(0,1),"outplot")
  hybrid <- mapa.choice(hybrid,c(TRUE,FALSE),"hybrid")
  paral <- mapa.choice(paral,c(0,1,2),"paral")
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Find locations in mapafit
  # Convert mapafit from earlier versions
  mapafit <- mapafit.convert(mapafit)
  
  # Get settings from mapafit
  ALs <- mapafit$AL[mapafit$use]
  maximumAL <- max(ALs)
  ets.type <- mapafit$type
  ppy <- mapafit$ppy
  
  # Default horizon. Note that fh = 0 is a valid input, giving only 
  # in-sample forecasts, so it cannot be used to request the default.
  if (is.null(fh) || fh == -1){
    fh <- ppy
  }
  
  # Override ifh.c to be >= fh if conf.lvl is not NULL
  # Output will crop values to ifh
  ifh.c <- ifh
  if (!is.null(conf.lvl) && fh > 0){
    if (ifh.c < fh){
      ifh.c <- fh
    } 
  } else {
    # Prediction intervals need out-of-sample forecasts
    conf.lvl <- NULL
  }
  
  # In-sample MAPA
  if (ets.type == "ets"){
    # ETS works with a single observation
    i.start <- max(ppy,maximumAL)
  } else if (ets.type == "es"){
    # ES needs at least 2 observations
    i.start <- max(ppy,maximumAL*2)
  }
  if (ifh.c>0 && i.start<observations){   # Do not produce in-sample forecasts if there is 
    # not enough sample
    infor <- array(NA,c(ifh.c,observations),dimnames=list(paste("t+",1:ifh.c,sep="")))
    if (paral != 0){  # Parallel run
      clout <- mapa.cluster(paral,cl)
      if (clout$created){
        on.exit(stopCluster(clout$cl))
      }
      inforlist <- clusterApplyLB(clout$cl, i.start:(observations-1), mapafor.loop, 
                                  y=y, mapafit=mapafit, ifh.c=ifh.c, comb=comb, 
                                  hybrid=hybrid, xreg=xreg)
      infor[, (i.start+1):observations] <- do.call(cbind,inforlist)
    } else {          # Serial run
      for (i in i.start:(observations-1)){
        infor[, i+1] <- mapafor.loop(i, y, mapafit, ifh.c, comb, hybrid, xreg)
      }    
    }
  } else {
    infor <- NULL
    ifh.c <- 0        # Override in-sample output if no forecasts are calculated
    conf.lvl <- NULL  # Confidence intervals cannot be calculated
    
  }
  
  # Out-of-sample MAPA
  if (fh>0){
    outfor <- mapacalc(y, mapafit, fh=fh, comb=comb, outplot=0, hybrid=hybrid, xreg=xreg)$forecast
  } else {
    outfor <- NULL
  }
  
  # Calculate in-sample errors
  if (ifh.c == 1) {
    resid <- y - t(infor)
    MSE <- array(mean(resid^2, na.rm=TRUE),c(1,1),dimnames=list("t+1","MSE"))
    MAE <- array(mean(abs(resid), na.rm=TRUE),c(1,1),dimnames=list("t+1","MAE"))
  } else if (ifh.c > 1) {
    MSE <- array(NA,c(ifh.c,1),dimnames=list(paste("t+",1:ifh.c,sep=""),"MSE"))
    MAE <- array(NA,c(ifh.c,1),dimnames=list(paste("t+",1:ifh.c,sep=""),"MAE"))
    for (h in 1:min(ifh.c,(observations-ppy))) {
      resid <- y[h:observations] - infor[h, 1:(observations-h+1)]
      MSE[h] <- mean(resid^2, na.rm=TRUE)
      MAE[h] <- mean(abs(resid), na.rm=TRUE)
    }
  } else {
    MSE <- NULL
    MAE <- NULL
  }
  
  # Calculate prediction intervals
  if (!is.null(conf.lvl)){
    intv.idx <- !is.na(MSE)
    if (sum(!intv.idx)>0){
      intv <- c(sqrt(MSE[intv.idx]), sqrt(MSE[1])*sqrt((sum(intv.idx)+1):fh))
    } else {
      intv <- c(sqrt(MSE[1:fh]))
    }
    conf.lvl[conf.lvl>0.999999999999999] <- 0.999999999999999
    conf.lvl <- unique(conf.lvl)
    PIn <- length(conf.lvl)
    conf.lvl <- sort(conf.lvl,decreasing=TRUE)
    z <- abs(qnorm((1-conf.lvl)/2))
    PI <- array(NA,c(2*PIn,fh))
    for (i in 1:PIn){
      PI[i,] <- outfor + intv*z[i]
      PI[2*PIn-i+1,] <- outfor - intv*z[i]
    }
    rownames(PI) <- c(paste("Upper",format(conf.lvl,digit=2)),
                      paste("Lower",format(conf.lvl[PIn:1],digit=2)))
    colnames(PI) <- paste("t+",1:fh,sep="")
  } else {
    PI <- NULL
  }
  
  # Crop insample forecasts and erros to ifh
  if (ifh.c > 0 && ifh > 0){
    infor <- array(infor[1:ifh,],c(ifh,observations),dimnames=list(paste("t+",1:ifh,sep="")))
    MSE <- MSE[1:ifh,]
    MAE <- MAE[1:ifh,]
  } else {
    infor <- NULL
    MSE <- NULL
    MAE <- NULL
  }
  
  # Construct output
  out <- structure(list("infor"=infor,"outfor"=outfor,"PI"=PI,"MSE"=MSE,"MAE"=MAE,
                        "y"=y,"fh"=fh,"ifh"=ifh,"comb"=comb,"hybrid"=hybrid,
                        "conf.lvl"=conf.lvl,"AL"=ALs,"ppy"=ppy,"type"=ets.type),class="mapa.for")
  
  # Produce plot
  if (outplot==1){
    plot(out)
  }
  
  return(out)
  
}

#-------------------------------------------------
mapafor.loop <- function(i, y, mapafit, ifh.c, comb, hybrid, xreg){
  # Internal function for producing the in-sample forecasts from origin i
  
  observations <- length(y)
  inobs <- as.matrix(y[1:i])
  infor <- mapacalc(inobs, mapafit, fh=ifh.c, comb=comb, outplot=0, hybrid=hybrid, xreg=xreg)$forecast
  # Crop out-of-sample predictions
  if ((i+ifh.c)>observations){
    k <- (i+ifh.c) - observations
    infor[(ifh.c-k+1):ifh.c] <- rep(NA,k)
  }
  
  return(infor)
  
}

#-------------------------------------------------
#' @rdname summary.mapa.fit
#' @export
print.mapa.for <- function(x,...){
  
  cat(paste0("MAPA forecasts (",x$type,", ",x$comb,ifelse(x$hybrid,", hybrid",""),")\n"))
  if (!is.null(x$MSE)){
    cat(paste0("In-sample t+1 MSE: ",round(x$MSE[1],2), ", MAE: ", round(x$MAE[1],2),"\n"))
  }
  if (is.null(x$outfor)){
    cat("No out-of-sample forecasts.\n")
  } else {
    cat("Out-of-sample forecasts:\n")
    if (is.null(x$PI)){
      fc <- x$outfor
      if (!isa(fc,"ts")){
        names(fc) <- paste("t+",1:length(fc),sep="")
      }
      print(fc)
    } else {
      # Order intervals from the widest lower to the widest upper
      sz <- dim(x$PI)[1]
      temp <- t(rbind(x$PI[sz:(1+sz/2),,drop=FALSE],as.numeric(x$outfor),x$PI[(sz/2):1,,drop=FALSE]))
      colnames(temp)[sz/2 + 1] <- "Forecast"
      print(temp)
    }
  }
  
  return(invisible(x))
  
}

#-------------------------------------------------
#' @rdname summary.mapa.fit
#' @export
summary.mapa.for <- function(object,...){
  # In-sample errors per horizon
  if (!is.null(object$infor)){
    y <- as.numeric(object$y)
    n <- length(y)
    ifh <- dim(object$infor)[1]
    err <- array(NA,c(ifh,4),dimnames=list(paste("t+",1:ifh,sep=""),
                                          c("N","ME","MAE","RMSE")))
    for (h in 1:ifh){
      resid <- y[h:n] - object$infor[h, 1:(n-h+1)]
      err[h,"N"] <- sum(!is.na(resid))
      err[h,"ME"] <- mean(resid, na.rm=TRUE)
      err[h,"MAE"] <- mean(abs(resid), na.rm=TRUE)
      err[h,"RMSE"] <- sqrt(mean(resid^2, na.rm=TRUE))
    }
  } else {
    err <- NULL
  }
  
  # Out-of-sample forecasts and intervals
  if (!is.null(object$outfor)){
    fc <- cbind("Forecast"=as.numeric(object$outfor))
    if (!is.null(object$PI)){
      sz <- dim(object$PI)[1]
      fc <- t(rbind(object$PI[sz:(1+sz/2),,drop=FALSE],"Forecast"=as.numeric(object$outfor),
                    object$PI[(sz/2):1,,drop=FALSE]))
    }
    rownames(fc) <- paste("t+",1:object$fh,sep="")
  } else {
    fc <- NULL
  }
  
  out <- structure(list("errors"=err,"forecasts"=fc,"comb"=object$comb,"hybrid"=object$hybrid,
                        "AL"=object$AL,"type"=object$type,"ppy"=object$ppy,"fh"=object$fh,
                        "conf.lvl"=object$conf.lvl,"observations"=length(object$y)),
                   class="summary.mapa.for")
  
  return(out)
  
}

#-------------------------------------------------
#' @export
print.summary.mapa.for <- function(x,digits=max(3,getOption("digits")-3),...){
  
  cat(paste0("MAPA forecasts using ",x$type,"\n"))
  cat(paste0("Observations: ",x$observations,"\tOriginal frequency: ",x$ppy,"\n"))
  cat(paste0("Aggregation levels: ",paste(x$AL,collapse=", "),"\n"))
  cat(paste0("Combination: ",x$comb,"\tHybrid: ",x$hybrid,"\tHorizon: ",x$fh,"\n"))
  if (!is.null(x$conf.lvl)){
    cat(paste0("Prediction intervals: ",paste(paste0(100*sort(x$conf.lvl),"%"),collapse=", "),
               " (empirical, from in-sample errors)\n"))
  }
  if (!is.null(x$errors)){
    cat("\nIn-sample errors by horizon:\n")
    print(round(x$errors,digits))
  }
  if (!is.null(x$forecasts)){
    cat("\nOut-of-sample forecasts:\n")
    print(round(x$forecasts,digits))
  }
  
  return(invisible(x))
  
}

#-------------------------------------------------
#' Plot MAPA in- and out-of-sample forecasts
#'
#' @description The following function plots the in- and out-of-sample forecasts and prediction intervals produced by mapafor or mapa. The plot is the same as the one that can be produced optionally by mapafor.
#'
#' @param x MAPA forecasts (from mapafor or mapa).
#' @param insample Plot the in-sample forecasts, if available.
#' @param PI Plot the prediction intervals, if available.
#' @param col Vector of 2 colours for the time series and the forecasts. If NULL the default colours are used.
#' @param col.PI Colour for the prediction intervals. When multiple intervals are available these are shaded from this colour. If NULL, a lighter shade of the forecast colour is used, so the intervals follow any colour set with \code{col}.
#' @param ... Additional arguments passed to \code{plot}, for example \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim}, and \code{lwd}.
#'
#' @note If y is a ts object, its time index is used for the horizontal axis. The current graphical layout is respected, so the plot can be placed in a user defined grid.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @seealso \code{\link{mapafor}}, \code{\link{mapa}}.
#'
#' @keywords ~plot
#' @examples
#' out <- mapa(admissions, ifh=3, conf.lvl=c(0.8,0.95))
#' plot(out)
#' plot(out, insample=FALSE, main="Admissions", ylab="Patients", col.PI="lightblue")
#' plot(out, PI=FALSE, lwd=2, col=c("grey50","blue"))
#'
#' # The prediction intervals are shaded from the forecast colour, unless col.PI is given
#' plot(out, col=c("grey30","darkgreen"))
#' @export
plot.mapa.for <- function(x,insample=c(TRUE,FALSE),PI=c(TRUE,FALSE),col=NULL,col.PI=NULL,...){
  insample <- insample[1]
  PI <- PI[1]
  
  y <- x$y
  observations <- length(y)
  fh <- x$fh
  outfor <- if (is.null(x$outfor)) NULL else as.numeric(x$outfor)
  infor <- if (insample) x$infor else NULL
  PIs <- if (PI) x$PI else NULL
  ifh <- if (is.null(infor)) 0 else dim(infor)[1]
  
  # Colours
  if (is.null(col)){
    col <- brewer.pal(3,"Set1")[c(2,1)]
  }
  col <- rep(col,length.out=2)
  if (is.null(col.PI)){
    # Shade the prediction intervals from the forecast colour, lightened towards white
    # col.PI <- rgb(1,0.75,0.75,1)
    col.PI <- 1-0.25*(1-c(col2rgb(col[2])/255))
    col.PI <- rgb(col.PI[1],col.PI[2],col.PI[3])
  }
  
  # Time index, extended to cover in-sample forecasts that end out-of-sample
  tm <- mapa.time(y,max(fh,ifh))
  tall <- c(tm$tin,tm$tout)
  tout <- tm$tout[seq_len(fh)]
  
  # Find min max
  args <- list(x=tm$tin, y=as.numeric(y), type="l", col=col[1], xlab="", ylab="", 
               main="Forecast", xlim=range(c(tm$tin,tout)), 
               ylim=mapa.ylim(y,outfor,PIs))
  args <- modifyList(args,list(...))
  do.call(plot,args)
  lwd <- if (is.null(args$lwd)) 1 else args$lwd
  
  # In-sample
  if (ifh>0){
    if (ifh==1){
      lines(tm$tin,infor[1,],col=col[2],lwd=lwd)
    } else {
      # clrs = rainbow(observations-ppy)
      for (i in which(colSums(!is.na(infor))>0)){
        lines(tall[i:(i+ifh-1)],infor[,i],col=col[2],lwd=lwd)
      }
    }
  }
  
  # Prediction intervals
  if (!is.null(PIs)){
    PIn <- dim(PIs)[1]/2
    for (i in 1:PIn){
      rr <- 1-(1-(i-1)*(0.5/PIn))*(1-c(col2rgb(col.PI)/255))
      polygon(c(tout,rev(tout)),c(PIs[i,],PIs[PIn*2+1-i,fh:1]),
              col=rgb(rr[1],rr[2],rr[3]), border=NA)
    }
  }

  # Out-of-sample
  if (!is.null(outfor)){
    if (ifh == 0){
      lines(c(tm$tin[observations],tout),c(y[observations],outfor),col=col[2],lwd=lwd)
    } else if (ifh == 1){
      lines(c(tm$tin[observations],tout),c(infor[1,observations],outfor),col=col[2],lwd=lwd)
    } else {
      lines(tout,outfor,col=col[2],lwd=lwd)
    }
  }
  
  return(invisible())
  
}
