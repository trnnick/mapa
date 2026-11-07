# Now split into multiple files (mapa.R, mapa.estim.R, mapa.forecast.R, mapa.utils.R)

#-------------------------------------------------
#' Multiple Aggregation Prediction Algorithm (Wrapper)
#'
#' @description
#' Wrapper to estimate and produce MAPA in- and out-of-sample forecasts.
#' Uses mapaest and mapafor.
#'
#' @param y In sample observations of a time series (vector). If y == "paper" then it prints paper reference.
#' @param ppy Periods in a season of the time series at the sampled frequency.
#'   If insample is a ts object then this is taken from its frequency, unless overriden.
#' @param fh Forecast horizon. Default = ppy. If fh = 0 then only in-sample forecasts are produced.
#' @param ifh In-sample forecast horizon. Default = 1.
#' @param minimumAL Lowest aggregation level to use. Default = 1.
#' @param maximumAL Highest aggregation level to use. Default = ppy, maximumAL>1.
#' @param comb Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
#' @param paral Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes and initialise cluster. Default is 0. Parallel processing is used for both the estimation and the in-sample forecasts.
#' @param display Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#' @param outplot Provide output plot. 0 = no; 1 = yes. Default is 0.
#' @param hybrid Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
#' @param model Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M", or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X", or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X", or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive, and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level.
#' @param type What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
#' @param conf.lvl Vector of confidence level for prediction intervals. Values must be (0,1). If conf.lvl == NULL then no intervals are calculated. For example to get the intervals for 80\% and 95\% use conf.lvl=c(0.8,0.95).
#' @param xreg Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample. Additional observations are unused. Note that including xreg will force type="es".
#' @param pr.comp MAPAx can use principal component analysis to preprocess xreg. When comp is -1 then the number of retained components is chosen automatically. When comp=0 then no pre-processing is performed and the original xreg is used. Any other value represents the number of principal components retained.
#' @param cl Cluster to use when paral=1, created with \code{parallel::makeCluster}. If NULL, a cluster named \code{cl} in the global environment is used, for compatibility with earlier versions.
#' @param ... Pass additional arguments to es or ets.
#'
#' @return
#' \item{infor}{In-sample forecasts.}
#' \item{outfor}{Out-of-sample forecasts.}
#' \item{PI}{Prediction intervals for given confidence levels.}
#' \item{MSE}{In-sample MSE error.}
#' \item{MAE}{In-sample MAE error.}
#' The output is of class \code{mapa.for}, which also stores the time series and the settings used. See \code{\link{mapafor}} and \code{\link{plot.mapa.for}}.
#'
#' @details This function calls \emph{mapaest} and \emph{mapafor} internally.
#'
#' Use \code{type="es"} for the \code{es} function of the \pkg{smooth} package, instead of the default \code{ets} function of the \pkg{forecast} package. See \code{\link{mapaest}} for details on the differences between the two and on the use of xreg.
#'
#' @note The calculation of the prediction intervals is based on the empirical multiple step ahead MSE. To speed up calculations set \code{conf.lvl=NULL}. If very long forecast horizons are requested then once no more t+h MSE can be calculated the following approximation is used: \eqn{\sqrt{MSE_{t+1}}*\sqrt{h}}{sqrt(MSE_[t+1])*sqrt(h)} for the error.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}; Fotios Petropoulos.
#'
#' @seealso \code{\link{mapaest}}, \code{\link{mapafor}}, \code{\link{mapasimple}}.
#'
#' @keywords ~mapaest ~mapafor ~mapasimple
#' @examples
#' # MAPA with ets from the forecast package (default), forecasting 18 steps 
#' # ahead, with in-sample forecasts up to 3 steps ahead and prediction intervals
#' out <- mapa(admissions, fh=18, ifh=3, conf.lvl=c(0.8,0.95))
#' print(out)
#' summary(out)
#' plot(out)
#'
#' # MAPAx with es from the smooth package. xreg must contain the in-sample 
#' # observations and the fh future values
#' set.seed(1)
#' fh <- 12
#' x <- matrix(rnorm((length(admissions)+fh)*2), ncol=2)
#' out <- mapa(admissions, fh=fh, type="es", xreg=x)
#'
#' \donttest{
#' # es with a different combination operator and without hybrid forecasts
#' out <- mapa(admissions, type="es", comb="median", hybrid=FALSE)
#'
#' # Non-seasonal models only, at aggregation levels 1 to 6
#' out <- mapa(admissions, model="ZZN", maximumAL=6)
#'
#' # MAPAx with xreg preprocessed by principal components, selecting the number
#' # of components automatically, and plotting the output
#' out <- mapa(admissions, fh=fh, type="es", xreg=x, pr.comp=-1, conf.lvl=0.9, outplot=1)
#'
#' # Only in-sample forecasts
#' out <- mapa(admissions, fh=0, ifh=6)
#' }
#' @export
mapa <- function(y, ppy=NULL, fh=NULL, ifh=1, minimumAL=1, maximumAL=ppy, 
	comb=c("w.mean","w.median","mean","median","wght"), paral=c(0,1,2), display=c(0,1), outplot=c(0,1), 
	hybrid=c(TRUE,FALSE), model="ZZZ", type=c("ets","es"), conf.lvl=NULL, 
	xreg=NULL, pr.comp=0, cl=NULL, ...){
  # Paper info
  if (!is.numeric(y)){
    return(mapa.paper())
  }
  
  # Defaults
  comb <- match.arg(comb,c("w.mean","w.median","mean","median","wght"))
  paral <- mapa.choice(paral,c(0,1,2),"paral")
  display <- mapa.choice(display,c(0,1),"display")
  outplot <- mapa.choice(outplot,c(0,1),"outplot")
  hybrid <- mapa.choice(hybrid,c(TRUE,FALSE),"hybrid")
  type <- match.arg(type,c("ets","es"))
  
  # Setup parallel processing if required
  # A single cluster is used for both estimation and forecasting
  if (paral == 2){
    clout <- mapa.cluster(paral,cl)
    cl <- clout$cl
    paral <- 1
    on.exit(stopCluster(cl))
  }
  
  # Estimate MAPA
  mapafit <- mapaest(y=y, ppy=ppy, minimumAL=minimumAL, 
                     maximumAL=maximumAL, paral=paral, display=display, 
                     outplot=outplot, model=model, type=type, xreg=xreg, pr.comp=pr.comp, 
                     cl=cl, ...)
  
  # If not horizon was given forecasts a full season
  # This is handled by mapafor, when fh=NULL
  
  # Produce in- and out-of-sample forecasts
  out <- mapafor(y=y, mapafit=mapafit, fh=fh, ifh=ifh, comb=comb, 
                 outplot=outplot, hybrid=hybrid, conf.lvl=conf.lvl, xreg=xreg,
                 paral=paral, cl=cl)
  
  return(out)
  
}

#-------------------------------------------------
mapacomb <- function(ALs,ppy,FCs,comb){
# This function combines the translated ets states
#
# Inputs:
#   ALs         = Aggregation levels that are combined. These correspond to the rows of FCs.
#   ppy         = Periods in a season at the original sampling frequency.
#   FCs         = Array of translated forecasts (aggregation levels x components x horizon).
#   comb        = Combination operator.
#
# Output:
#   forecasts   = Combined forecasts.
#   combined    = Combined components (components x horizon).
#   perm_seas   = Aggregation levels for which seasonality is permitted.
  
  # perm_levels is not needed for forecasting. This is already checked in the estimation.
  # Aggregation levels that were not estimated are not in FCs.
  # perm_levels <- array(0, maximumAL) # permitted levels due to ETS implementation (observations>=4)
  # if (observations %/% AL >=4){
  #   perm_levels[AL] <- 1
  # }
  # perm_levels <- perm_levels[minimumAL:maximumAL]
  perm_seas <- (ppy %% ALs == 0) & (ALs<ppy)   # permitted seasonalities
  fh <- dim(FCs)[3]
  
  # Combination operator for level, trend, season, and xreg
  operator <- switch(comb,
                     "mean"=rep("mean",4), # alternative averaging operators
                     "median"=rep("median",4),
                     "wght"=rep("wght",4), # Weighted sum
                     "w.mean"=c("mean","mean","wght","wght"),
                     "w.median"=c("median","median","wght","wght"))
  
  # Combine each component
  # Forecast multiple steps ahead, or a single step ahead, with the same calculations
  # Seasonality is combined only across the levels that it is permitted
  cmp.names <- c("Level","Trend","Season","Xreg")
  combined <- array(0,c(4,fh),dimnames=list(cmp.names,dimnames(FCs)[[3]]))
  for (k in 1:4){
    if (cmp.names[k] == "Season"){
      idx <- perm_seas
    } else {
      idx <- rep(TRUE,length(ALs))
    }
    # Check that all are arrays
    X <- matrix(FCs[idx, k+1, ],nrow=sum(idx),ncol=fh)
    combined[k,] <- mapacomb.comp(X,ALs[idx],operator[k])
  }
  
  # MAPA forecasts
  forecasts <- colSums(combined,na.rm=TRUE)
  
  # Return output
  return(list(forecasts=forecasts,combined=combined,perm_seas=perm_seas))
}

#-------------------------------------------------
mapacomb.comp <- function(X,ALs,operator){
# Combine a single component across aggregation levels
#
# Inputs:
#   X           = Matrix of a component (aggregation levels x horizon).
#   ALs         = Aggregation levels corresponding to the rows of X.
#   operator    = "mean", "median", or "wght" - where each aggregation level 
#                 is weighted inversly to aggregation. 
  
  # No levels to combine, so the component does not contribute
  if (dim(X)[1] == 0){
    return(rep(0,dim(X)[2]))
  }
  
  if (operator == "mean"){
    out <- colMeans(X) # MAPA(mean) forecasts
  } else if (operator == "median"){
    # forecasts <- colSums(rbind(colMedians(level),colMedians(trend),
    #                            colMedians(season)),na.rm=TRUE) # MAPA(median) forecasts
    out <- apply(X,2,"median") # MAPA(median) forecasts
  } else if (operator == "wght"){
    # Weighted sum
    wghts <- 1/ALs
    wghts <- wghts/sum(wghts)
    out <- colSums(X * wghts)
  }
  
  return(out)
  
}

#-------------------------------------------------
#' MAPA Single Origin Forecast (Wrapper)
#'
#' @description
#' The following function estimates MAPA and produces out-of-sample forecasts.
#' It is not suggested, unless only simple single origin forecasts are required.
#'
#' @param y    In sample observations of a time series (vector). If y == "paper" then
#'   it prints paper reference.
#' @param ppy Periods in a season of the time series at the sampled frequency.
#'   If insample is a ts object then this is taken from its frequency, unless overriden.
#' @param fh Forecast horizon. Default = ppy.
#' @param minimumAL Lowest aggregation level to use. Default = 1, maximumAL>1.
#' @param maximumAL Highest aggregation level to use. Default = ppy.
#' @param comb Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
#' @param paral    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes
#'   and initialise cluster. Default is 0.
#' @param display Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#' @param outplot Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, forecasts, and components. For the components the spectral colouring scheme is used. Red is the lowest aggregation level and blue the highest. Default is 0.
#' @param hybrid Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
#' @param model Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M", or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X", or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X", or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive, and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. For aggregation levels that seasonality becomes 1 then a non-seasonal model is estimated.
#' @param type What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
#' @param xreg Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample. Additional observations are unused. Note that including xreg will force type="es".
#' @param pr.comp MAPAx can use principal component analysis to preprocess xreg. When comp is -1 then the number of retained components is chosen automatically. When comp=0 then no pre-processing is performed and the original xreg is used. Any other value represents the number of principal components retained.
#' @param cl Cluster to use when paral=1, created with \code{parallel::makeCluster}. If NULL, a cluster named \code{cl} in the global environment is used, for compatibility with earlier versions.
#' @param ... Pass additional arguments to es or ets.
#'
#' @return
#' \item{forecast}{Vector with forecasts. If y is a ts object, then this is also a ts object.}
#' \item{components}{Array with MAPA components.}
#' \item{combined}{Combined components across aggregation levels.}
#' The output is of class \code{mapa.calc}. See \code{\link{mapacalc}} and \code{\link{plot.mapa.calc}}.
#'
#' @details This function calls \emph{mapaest} and \emph{mapacalc} internally.
#'
#' Use \code{type="es"} for the \code{es} function of the \pkg{smooth} package, instead of the default \code{ets} function of the \pkg{forecast} package. See \code{\link{mapaest}} for details on the differences between the two and on the use of xreg.
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
#' @seealso \code{\link{mapa}}.
#'
#' @keywords ~mapa
#' @examples
#' # MAPA forecasts with ets from the forecast package (default)
#' out <- mapasimple(admissions)
#' print(out)
#'
#' # With es from the smooth package, 24 steps ahead, plotting forecasts and components
#' out <- mapasimple(admissions, fh=24, type="es", outplot=2)
#'
#' # Weighted combination, without the hybrid forecasts
#' out <- mapasimple(admissions, comb="wght", hybrid=FALSE)
#'
#' # MAPAx: xreg must contain the in-sample observations and the fh future values
#' set.seed(1)
#' x <- rnorm(length(admissions)+12)
#' out <- mapasimple(admissions, fh=12, type="es", xreg=x)
#' summary(out)
#' @export
mapasimple <- function(y, ppy=NULL, fh=NULL, minimumAL=1, maximumAL=ppy, comb=c("w.mean","w.median","mean","median","wght"), 
                       paral=c(0,1,2), display=c(0,1), outplot=c(0,1,2), hybrid=c(TRUE,FALSE), 
                       model="ZZZ", type=c("ets","es"), xreg=NULL, pr.comp=0, cl=NULL, ...){
  # Paper info
  if (!is.numeric(y)){
    return(mapa.paper())
  }  
  
  # Defaults
  comb <- match.arg(comb,c("w.mean","w.median","mean","median","wght"))
  paral <- mapa.choice(paral,c(0,1,2),"paral")
  display <- mapa.choice(display,c(0,1),"display")
  outplot <- mapa.choice(outplot,c(0,1,2),"outplot")
  hybrid <- mapa.choice(hybrid,c(TRUE,FALSE),"hybrid")
  type <- match.arg(type,c("ets","es"))
  
  # Estimate MAPA
  mapafit <- mapaest(y=y, ppy=ppy, minimumAL=minimumAL, 
                     maximumAL=maximumAL, paral=paral, display=display, 
                     outplot=outplot, model=model, type=type, xreg=xreg, pr.comp=pr.comp, 
                     cl=cl, ...)
  
  # If not horizon was given forecasts a full season
  # This is handled by mapacalc, when fh=NULL
  
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
  #   pr.comp    = list containing list("pr.comp", "mean", "sd", "rotation"). If 
  #                rotation is provided, as is done when forecasting, it is reused 
  #                instead of calculating the principal components again. This keeps 
  #                the components consistent with the ones used in estimation.
  
  p <- dim(x)[2]
  
  # Get bits from list
  comp <- pr.comp$pr.comp
  x.mn <- pr.comp$mean
  x.sd <- pr.comp$sd
  x.rot <- pr.comp$rotation
  
  if (comp != 0){
    # Reduction is perfomed
    
    # Normalise xreg
    n <- dim(x)[1]
    x <- (x - matrix(rep(x.mn,n),nrow=n,byrow=TRUE))/matrix(rep(x.sd,n),nrow=n,byrow=TRUE)
    
    # Check if requested number of components is possible
    if (comp > p){
      stop("The requested number of principal components exceeds the number of xreg variables.")
    }
    
    if (!is.null(x.rot)){
      # Use the rotation from estimation
      x.out <- x %*% x.rot
    } else {
    
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
      x.rot <- x.pr$rotation[,1:comp,drop=FALSE]
      
    }
    
  } else {
    # No pre-processing
    x.out <- x
  }
  
  return(list("x.out"=x.out,"pr.comp"=list("pr.comp"=comp,"mean"=x.mn,"sd"=x.sd,"rotation"=x.rot)))
  
}