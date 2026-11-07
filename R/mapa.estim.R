# These functions are handling the estimation of MAPA
# Nikolaos Kourentzes 2016 

#' MAPA Estimation
#'
#' @description The following function estimates MAPA and saves all fitted model.
#'
#' @param y In sample observations of a time series (vector).
#' @param ppy Periods in a season of the time series at the sampled frequency.
#'   If insample is a ts object then this is taken from its frequency, unless overriden.
#' @param minimumAL Lowest aggregation level to use. Default = 1, maximumAL>1.
#' @param maximumAL Highest aggregation level to use. Default = ppy.
#' @param paral    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes
#'   and initialise cluster. Default is 0.
#' @param display Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#' @param outplot Provide output plot. 0 = no; 1 or 2 = yes. Default is 0.
#' @param model Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M", or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X", or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X", or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive, and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. For aggregation levels that seasonality becomes 1 then a non-seasonal model is estimated.
#' @param type What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
#' @param xreg Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample. Additional observations are unused. Note that including xreg will force type="es".
#' @param pr.comp MAPAx can use principal component analysis to preprocess xreg. When comp is -1 then the number of retained components is chosen automatically. When comp=0 then no pre-processing is performed and the original xreg is used. Any other value represents the number of principal components retained.
#' @param cl Cluster to use when paral=1, created with \code{parallel::makeCluster}. If NULL, a cluster named \code{cl} in the global environment is used, for compatibility with earlier versions.
#' @param ... Pass additional arguments to es or ets.
#'
#' @return
#' Estimated MAPA model structure, of class \code{mapa.fit}. This is a list containing:
#' \item{models}{List with the fitted model at each aggregation level, as returned by \code{ets} or \code{es}. NULL for levels that could not be estimated.}
#' \item{AL}{Aggregation levels.}
#' \item{use}{TRUE for the aggregation levels that were estimated.}
#' \item{q}{Sample size at each aggregation level.}
#' \item{ppyA}{Seasonal periods at each aggregation level.}
#' \item{ppy, type, model, observations}{Settings used in the estimation.}
#' \item{xreg.n}{Number of xreg variables.}
#' \item{pr.comp}{Information for the preprocessing of xreg.}
#' Earlier versions (<= 2.0.6) stored the fitted models in a list-matrix. Such objects are still accepted by \code{mapacalc}, \code{mapafor}, and the methods for \code{mapa.fit}.
#'
#' @details
#' MAPA can use two implementations of exponential smoothing at each aggregation level. With \code{type="ets"} the \code{ets} function of the \pkg{forecast} package is used, which supports seasonal periods up to 24. With \code{type="es"} the \code{es} function of the \pkg{smooth} package is used, which supports any seasonal period, the model codes "X" and "Y", and exogenous variables (xreg). Additional arguments are passed to \code{ets} or \code{es} through \code{...}, for example the information criterion used for model selection (\code{ic="bic"} for \code{ets}, or \code{ic="BIC"} for \code{es}).
#'
#' MAPAx: exogenous variables are provided with \code{xreg}, as a vector or as a matrix with one column per variable. These are aggregated in the same way as y at each aggregation level. Including xreg forces \code{type="es"}. For estimation xreg must have at least as many observations as y, and any additional observations are ignored. For forecasting, with \code{mapacalc}, \code{mapafor}, \code{mapa}, and \code{mapasimple}, xreg must also contain the values for the forecast horizon, so it must have at least length(y) + fh observations. The same xreg variables must be used in estimation and forecasting. When \code{pr.comp} is not 0, xreg is preprocessed with principal components. These are calculated in estimation and stored in the fitted MAPA, so that the same components are used in forecasting.
#'
#' Since version 2.1.0 the fitted MAPA is a list, as described in Value. Earlier versions (<= 2.0.6) returned a list-matrix with one row per aggregation level. Code that depends on the earlier structure, for example indexing as \code{mapafit[,"use"]} or \code{mapafit[[1,14]]}, can transform the fitted MAPA to the earlier structure with \code{\link{mapafit.old}}. This is provided only for compatibility and will be removed in a later version of the package. Such code should be updated to use the current structure, for example \code{mapafit$use} or \code{mapafit$models[[1]]$components}.
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
#' @seealso \code{\link{mapa}}, \code{\link{mapafor}}.
#'
#' @keywords ~mapa ~mapafor
#' @examples
#' # Estimate MAPA using ets from the forecast package (default)
#' mapafit <- mapaest(admissions)
#' print(mapafit)
#' summary(mapafit)
#' plot(mapafit)
#'
#' # Use es from the smooth package instead
#' mapafit.es <- mapaest(admissions, type="es")
#'
#' # Restrict the model at each aggregation level. With type="es" the codes "X"
#' # and "Y" can be used, here to select between no and additive trend and 
#' # seasonality, with additive errors
#' mapafit.es <- mapaest(admissions, type="es", model="AXX")
#'
#' # Use only non-seasonal models, at aggregation levels 2 to 6
#' mapafit <- mapaest(admissions, model="ZZN", minimumAL=2, maximumAL=6)
#'
#' # If y is not a ts object, then ppy must be provided
#' mapafit <- mapaest(as.numeric(admissions), ppy=12)
#'
#' # MAPAx: estimate with exogenous variables, here 3 simulated ones. Only 
#' # type="es" accepts xreg. Observations beyond the length of y are ignored 
#' # in estimation, so the same xreg can be used later for forecasting
#' set.seed(1)
#' x <- matrix(rnorm((length(admissions)+12)*3), ncol=3)
#' mapafit.x <- mapaest(admissions, type="es", xreg=x)
#' summary(mapafit.x)
#' plot(mapafit.x)
#'
#' \donttest{
#' # Preprocess xreg with principal components, retaining 2 components, or
#' # selecting the number of components automatically
#' mapafit.pc <- mapaest(admissions, type="es", xreg=x, pr.comp=2)
#' mapafit.pc <- mapaest(admissions, type="es", xreg=x, pr.comp=-1)
#' summary(mapafit.pc)
#'
#' # Pass additional arguments to ets or es, here the information criterion
#' mapafit <- mapaest(admissions, ic="bic")
#' mapafit.es <- mapaest(admissions, type="es", ic="BIC")
#'
#' # Parallel processing with an existing cluster
#' library(parallel)
#' cl <- makeCluster(2)
#' mapafit <- mapaest(admissions, paral=1, cl=cl)
#' stopCluster(cl)
#' }
#'
#' # Code that depends on the structure of earlier versions can use mapafit.old.
#' # This will be removed in a later version of the package.
#' mapafit.prev <- suppressWarnings(mapafit.old(mapafit))
#' mapafit.prev[,"use"]
#' @export
mapaest <- function(y, ppy=NULL, minimumAL=1, maximumAL=ppy, paral=c(0,1,2), display=c(0,1), 
                    outplot=c(0,1,2), model="ZZZ", type=c("ets","es"), xreg=NULL, 
                    pr.comp=0, cl=NULL, ...) {
  # Default type
  paral <- mapa.choice(paral,c(0,1,2),"paral")
  display <- mapa.choice(display,c(0,1),"display")
  outplot <- mapa.choice(outplot,c(0,1,2),"outplot")
  type <- match.arg(type,c("ets","es"))
  
  # Get ppy and maximumAL
  if (is.null(ppy)){
    if (isa(y,"ts")){
      ppy <- frequency(y)
      if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
      stop(paste("Input ppy is not given and y input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }  
  
  # If xreg is used then force type="es" and give a warning
  if (!is.null(xreg) && type=="ets"){
    type <- "es"
    warning('Only type="es" can accept xreg inputs. Forcing appropriate type.')
  }
  
  # Convert xreg to array and check its size
  if (!is.null(xreg)){
    n <- length(y)
    # Force xreg to be a matrix, each column a variable
    xreg <- mapa.xreg(xreg)
    p <- dim(xreg)[2]
    # Check that it contains at least n observations and trim
    xn <- dim(xreg)[1]
    if (xn < n){
      stop("Number of observations of xreg must be equal of exceed the length of y.")
    } else {
      xreg <- xreg[1:n,,drop=FALSE]
    }
    
    # Perform preprocessing
    x.m <- colMeans(xreg)
    x.sd <- apply(xreg,2,sd)
    pr.comp <- list("pr.comp"=pr.comp, "mean"=x.m, "sd"=x.sd)
    pr.out <- mapaprcomp(xreg,pr.comp)
    xreg <- pr.out$x.out
    pr.comp <- pr.out$pr.comp
    
  } else {
    # No xreg, so mark no pre-processing
    p <- 0
    pr.comp <- list("pr.comp"=0, "mean"=NULL, "sd"=NULL)
  }
  
  # Suggest to use es instead of ets when ppy > 24
  if (ppy > 24 && type == "ets"){
    warning('Your data have potential seasonality longer than 24 periods. Consider switching to type="es".')
  }
  
  # Make sure that maximumAL > 1 and larger than minimumAL
  if (maximumAL == 1){
    maximumAL = maximumAL + 1
  }
  
  if (minimumAL>=maximumAL){
    stop("maximumAL must be larger than minimumAL")
  }
  
  # Setup parallel processing if required
  if (paral != 0){
    clout <- mapa.cluster(paral,cl)
    cl <- clout$cl
    if (clout$created){
      # Stop parallel processing, even if estimation fails
      on.exit(stopCluster(cl))
    }
  }
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Override model if used with ETS
  if (type=="ets"){
    mn <- nchar(model)
    model.set <- substring(model, seq(1,mn,1), seq(1,mn,1))
    idx <- ((model.set)=="X" | (model.set)=="Y")
    if (any(idx)){
      warning('When type=="ets" neither "X" or "Y" can be used in model specification. These will be changed to "Z".')
      model.set[idx] <- "Z"
      model <- paste0(model.set,collapse="")
    }
  }
  
  # Aggregation and estimation
  ALvec <- minimumAL:maximumAL
  if (paral != 0){  # Parallel run
    mapafit <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapaest.loop, 
                              y=y, minimumAL=minimumAL, maximumAL=maximumAL, observations=observations, ppy=ppy,
                              display=display,model=model,type=type,xreg=xreg,...)  
  } else {          # Serial run
    mapafit <- vector("list", (maximumAL-minimumAL+1))
    for (i in 1:(maximumAL-minimumAL+1)){
      mapafit[[i]] <- mapaest.loop(i, y, minimumAL, maximumAL, observations, 
                                   ppy, display,model=model,type=type,xreg=xreg,...)
    }
  }
  
  # Process output
  # Re-arrange output for clusterApplyLB function
  # Fitted models are kept complete, so that changes in the forecast and smooth 
  # packages do not affect the stored information.
  models <- lapply(mapafit,function(x){x$fit})
  names(models) <- paste("AL",ALvec,sep="")
  use <- unlist(lapply(mapafit,function(x){x$use}))
  if (!any(use)){
    stop("Sample is too short to estimate the requested model at any aggregation level.")
  }
  mapafit <- structure(list("models"=models, "AL"=ALvec, "use"=use,
                            "q"=unlist(lapply(mapafit,function(x){x$q})),
                            "ppyA"=unlist(lapply(mapafit,function(x){x$ppyA})),
                            "ppy"=ppy, "type"=type, "model"=model, "observations"=observations,
                            "xreg.n"=p, "pr.comp"=pr.comp), class = "mapa.fit")
  
  # Plot model selection summary
  if (outplot > 0){
    plot(mapafit)
  }
    
  # Return output
  return(mapafit)

}

#-------------------------------------------------
mapaest.loop <- function(ALi, y, minimumAL, maximumAL, observations, 
                         ppy, display, model, type, xreg=NULL, ...){ 
  # Internal function for running a single loop in mapaest
  
  # Create ETS model strings
  mn <- nchar(model)
  model.set <- substring(model, seq(1,mn,1), seq(1,mn,1))
  if (type == "ets"){
    if (model.set[2]=="Z"){
      ets.damped <- NULL
    } else if (mn==4){
      ets.damped <- TRUE
    } else {
      ets.damped <- FALSE
    }
  } else if (type == "es"){
    if (mn == 4){
      model.set <- c(model.set[1],paste0(model.set[2],model.set[3]),model.set[mn])
    }
  } else {
    stop("Unrecognised exponential smoothing implementation type.")
  }
  
  ALvec <- minimumAL:maximumAL
  AL <- ALvec[ALi]
  
  # Console output
  if (display==1){
    txtc <- paste("Aggregation level: ",AL,"/",maximumAL,
                  " (",round(100*ALi/(maximumAL-minimumAL+1),2),"%)",sep="")
    cat(txtc)
  }
  
  # Observations and periods per year in the aggregated level
  alinfo <- mapa.alinfo(AL,observations,ppy)
  q <- alinfo$q
  ppyA <- alinfo$ppyA
  qA <- q # Keep for output, as q may be overriden below
  
  # Check if requested ETS model is possible for current AL
  # Different treatment of damped trend and npars between packages
  if (type == "ets"){
    npars <- 2
    if (!is.null(ets.damped)){ 
      npars <- npars + as.numeric(ets.damped)}
  } else if (type == "es"){
    npars <- 3
    if (nchar(model.set[2]) == 2){ 
      npars <- npars + 3}
  }
  # Common amongst both types
  if (model.set[2] == "A" | model.set[2] == "M"){
    npars <- npars + 2}
  if (tail(model.set,1) == "A" | tail(model.set,1) == "M"){
    npars <- npars + ppyA}
  if (q <= npars){
    q <- 1} # This will not estimate current AL
  
  # Adjust for explanatory variables
  if (!is.null(xreg)){
    xreg.p <- dim(xreg)[2]
  } else {
    xreg.p <- 0
  }
  
  if (q >= 4 + xreg.p){ # Fit only when there is enough sample
    
    # Aggregation
    yA <- mapa.aggr(y,AL,q)
    ats <- ts(yA, frequency = ppyA) # create the time series
    
    # Aggregate xreg
    if (!is.null(xreg)){
      xregA <- mapa.aggr(xreg,AL,q)
    } else {
      xregA <- NULL
    }
    
    # Check if seasonality exists and select appropriate model
    if ((tail(model.set,1) == "A" | tail(model.set,1) == "M") && ppyA == 1){
      mapa.model <- paste0(model.set[1],model.set[2],"N")
    } else if ((ppyA == 2 && q <= 3*ppyA) && type == "ets"){   
      # This is added to deal with the bug in the fourier function that 
      # is used for the initialisation of ets
      mapa.model <- paste0(model.set[1],model.set[2],"N")
    } else {
      mapa.model <- paste0(model.set[1],model.set[2],tail(model.set,1))
    }
    
    # Fit exponential smoothing
    if (type == "ets"){
      
      fit <- ets(ats,model=mapa.model,damped=ets.damped,...)
      
      # If time series is constant then ets does not return the same
      # structure. As the complete fitted object is now stored, this 
      # no longer needs to be corrected. 
      
    } else if (type == "es"){
      
      # When xreg is a single column array es wants the input as a vector
      if (!is.null(xreg)){
        if (dim(xregA)[2] ==1){
          xregA <- as.vector(xregA)
        }
      }
      # Turn off warnings for es - this is done when the model reduces pool 
      # due to sample size.
      fit <- suppressWarnings(es(ats,model=mapa.model,silent=TRUE,xreg=xregA, ...))
            
    }
    
    # Common part
    use <- TRUE
    
  } else {
    # No fit is taking place
    fit <- NULL
    use <- FALSE
  }
  
  # Update console display
  if (display==1){
    nc = nchar(txtc)
    cat("\r",strrep(" ",nc),"\r",sep="")
  }
  
  # Return loop result
  return(list("fit"=fit, "use"=use, "AL"=AL, "q"=qA, "ppyA"=ppyA))
  
}

#-------------------------------------------------
#' Produce estimated MAPA fit plot (OBSOLETE!)
#'
#' @description
#' The following function plots the fitted components at each aggregation
#' level of an estimated MAPA. The plot is the same as the one that can be
#' produced optionally by mapaest.
#'
#' @param mapafit Fitted MAPA model (from mapaest).
#'
#' @note The function 'plotmapa' is deprecated. Use plot() instead.
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
#' @seealso \code{\link{mapaest}}, \code{\link{mapa}}.
#'
#' @keywords ~plotmapa
#' @examples
#' mapafit <- mapaest(admissions,outplot=0)
#' # plotmapa is deprecated, use plot instead
#' plot(mapafit)
#' @export
plotmapa <- function(mapafit){
  
  warning("'plotmapa' is deprecated.\n Use plot() instead.")
  plot.mapa.fit(mapafit)
  
}

#-------------------------------------------------
#' Transform MAPA fit to the structure of earlier versions
#'
#' @description
#' Since version 2.1.0 the fitted MAPA model returned by \code{mapaest} is a list that keeps the complete fitted model at each aggregation level. Earlier versions (<= 2.0.6) stored it as a list-matrix, with one row per aggregation level. This function transforms the current structure to the earlier one, so that code that depends on the old indexing of the fitted MAPA keeps working.
#'
#' This function is provided only for compatibility and will be removed in a later version of the package. A warning is given every time it is used.
#'
#' @param mapafit Fitted MAPA model (from mapaest).
#'
#' @return Fitted MAPA model as a list-matrix of class \code{mapa.fit}, with one row per aggregation level, as in MAPA <= 2.0.6. For \code{type="ets"} the columns are the 18 elements of the \code{ets} fit that were stored previously, followed by \code{use}, \code{AL}, \code{original.ppy}, \code{etstype}, and \code{pr.comp}. For \code{type="es"} the columns are the elements of the \code{es} fit that were stored previously, followed by \code{smoothLength}, \code{use}, \code{AL}, \code{original.ppy}, \code{etstype}, and \code{pr.comp}. The output is still accepted by \code{mapacalc}, \code{mapafor}, and the methods of \code{mapa.fit}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @seealso \code{\link{mapaest}}.
#'
#' @keywords ~mapaest
#' @examples
#' mapafit <- mapaest(admissions)
#'
#' # Transform to the structure of earlier versions. A warning is given, as this 
#' # function will be removed in a later version of the package.
#' mapafit.prev <- suppressWarnings(mapafit.old(mapafit))
#'
#' # Old style indexing
#' mapafit.prev[,"use"]
#' mapafit.prev[[1,"method"]]
#'
#' # The equivalent with the current structure
#' mapafit$use
#' mapafit$models[[1]]$method
#'
#' # The earlier structure is still accepted for forecasting
#' mapacalc(admissions, mapafit.prev)
#' @export
mapafit.old <- function(mapafit){
  warning(paste("mapafit.old() is provided only for compatibility with code written for",
                "MAPA <= 2.0.6 and will be removed in a later version of the package.",
                "Please update your code to use the current mapafit structure, see ?mapaest."),
          call.=FALSE)
  
  if (is.matrix(mapafit)){
    # Already in the old structure
    return(mapafit)
  }
  
  # Fields that were stored for each package
  if (mapafit$type == "ets"){
    fields <- c("loglik","aic","bic","aicc","mse","amse","fit","residuals","fitted",
                "states","par","m","method","components","call","initstate",
                "sigma2","x")
  } else if (mapafit$type == "es"){
    acceptres <- c('model', 'timeElapsed', 'data', 'holdout', 'fitted', 
                   'residuals', 'forecast', 'states', 'profile', 'profileInitial', 
                   'persistence', 'phi', 'transition', 'measurement', 'initial', 
                   'initialType', 'initialEstimated', 'orders', 'arma', 'constant', 
                   'nParam', 'occurrence', 'formula', 'regressors', 'loss', 
                   'lossValue', 'logLik', 'distribution', 'scale', 'other', 
                   'B', 'lags', 'lagsAll', 'FI', 'ICs', 'call', 'bounds')
    # Only the fields that exist in the fitted objects were kept, so that 
    # all rows have the same columns 
    fields <- acceptres[acceptres %in% names(mapafit$models[[which(mapafit$use)[1]]])]
  }
  
  # Construct one row per aggregation level
  k <- length(mapafit$AL)
  rows <- vector("list",k)
  for (i in 1:k){
    fit <- mapafit$models[[i]]
    # Missing fields, or levels that were not estimated, are NULL
    fit <- lapply(fields,function(nm){fit[[nm]]})
    names(fit) <- fields
    if (mapafit$type == "es"){
      if (mapafit$use[i]){
        fit$smoothLength <- length(fields)
      } else {
        fit["smoothLength"] <- list(NULL)
      }
    }
    fit$use <- mapafit$use[i]
    fit$AL <- mapafit$AL[i]
    fit$original.ppy <- mapafit$ppy
    fit$etstype <- mapafit$type
    fit$pr.comp <- mapafit$pr.comp
    rows[[i]] <- rbind(fit)
  }
  
  # Process output
  mapafit <- do.call(rbind, rows)
  rownames(mapafit) <- paste("AL",unlist(mapafit[,"AL"]),sep="")
  mapafit <- structure(mapafit, class = "mapa.fit")
  
  return(mapafit)
  
}

#-------------------------------------------------
#' @rdname summary.mapa.fit
#' @export
print.mapa.fit <- function(x,...){
  
  # Get basic settings
  x <- mapafit.convert(x)
  
  cat(paste0("MAPA fitted using ",x$type, "\tOriginal frequency: ",x$ppy,"\n"))
  for (i in 1:length(x$AL)){
    if (!x$use[i]){
      cat(paste0("Aggregation level: ",x$AL[i],"\t","Method: not estimated (insufficient sample)\n"))
      next
    }
    fit <- x$models[[i]]
    if (x$type == "ets"){
      mdl <- gsub(',','',fit$method)
      x.n <- 0
    } else if (x$type == "es"){
      mdl <- fit$model
      # Check for xreg
      x.n <- length(fit$initial$xreg)
    }
    if (x.n == 0){ # No xreg
      cat(paste0("Aggregation level: ",x$AL[i],"\t","Method: ",mdl,"\n"))
    }else {
      cat(paste0("Aggregation level: ",x$AL[i],"\t","Method: ",mdl,"+X(",x.n,")","\n"))
    }
  }
  
  return(invisible(x))
    
}

#-------------------------------------------------
#' Summaries of MAPA estimation and forecasts
#'
#' @description Print and summary methods for the outputs of mapaest, mapacalc, mapasimple, mapafor, and mapa.
#'
#' @param object,x Output of mapaest (class \code{mapa.fit}), mapacalc or mapasimple (class \code{mapa.calc}), or mapafor or mapa (class \code{mapa.for}).
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' The summary of \code{mapa.fit} reports for each aggregation level: the aggregated sample size, the seasonal periods, the fitted model, the smoothing parameters, and the number of xreg variables. Levels that could not be estimated due to insufficient sample are marked.
#'
#' The summary of \code{mapa.calc} reports the combined level, trend, season, and xreg components for each horizon, alongside the forecasts.
#'
#' The summary of \code{mapa.for} reports the in-sample errors (ME, MAE, and RMSE) for each in-sample forecast horizon, and the out-of-sample forecasts with any prediction intervals.
#'
#' The summaries are returned as objects, so their contents can be used directly, for example \code{summary(fit)$table}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @seealso \code{\link{mapaest}}, \code{\link{mapacalc}}, \code{\link{mapafor}}.
#'
#' @keywords ~summary
#' @examples
#' mapafit <- mapaest(admissions)
#' summary(mapafit)
#' summary(mapafit)$table
#'
#' # Summary of a MAPAx with principal component preprocessing of xreg
#' set.seed(1)
#' x <- matrix(rnorm(length(admissions)*3), ncol=3)
#' summary(mapaest(admissions, type="es", xreg=x, pr.comp=-1))
#'
#' # Summaries of forecasts
#' summary(mapacalc(admissions, mapafit))
#' out <- mapafor(admissions, mapafit, ifh=3, conf.lvl=0.9)
#' summary(out)
#' summary(out)$errors
#' @export
summary.mapa.fit <- function(object,...){
  object <- mapafit.convert(object)
  k <- length(object$AL)
  
  # Table with information per aggregation level
  tbl <- data.frame("AL"=object$AL, "n"=object$q, "m"=object$ppyA, 
                    "Model"=NA, "alpha"=NA, "beta"=NA, "gamma"=NA, "phi"=NA,
                    "Xreg"=0, stringsAsFactors=FALSE)
  
  for (i in 1:k){
    if (!object$use[i]){
      tbl$Model[i] <- "-"
      next
    }
    fit <- object$models[[i]]
    components <- mapa.components(fit,object$type)
    tbl$Model[i] <- paste0("ETS(",paste0(components,collapse=","),")")
    if (object$type == "ets"){
      prm <- fit$par
    } else if (object$type == "es"){
      prm <- fit$persistence
      if (substr(components[2],2,2) == "d"){
        prm <- c(prm,"phi"=fit$phi)
      }
      tbl$Xreg[i] <- length(fit$initial$xreg)
    }
    for (prmname in c("alpha","beta","gamma","phi")){
      if (prmname %in% names(prm)){
        tbl[i,prmname] <- prm[[prmname]]
      }
    }
  }
  
  out <- structure(list("table"=tbl, "type"=object$type, "ppy"=object$ppy, 
                        "model"=object$model, "observations"=object$observations,
                        "xreg.n"=object$xreg.n, "pr.comp"=object$pr.comp$pr.comp),
                   class="summary.mapa.fit")
  
  return(out)
  
}

#-------------------------------------------------
#' @export
print.summary.mapa.fit <- function(x,digits=3,...){
  
  cat(paste0("MAPA fitted using ",x$type,"\n"))
  if (!is.na(x$observations)){
    cat(paste0("Observations: ",x$observations,"\t"))
  }
  cat(paste0("Original frequency: ",x$ppy))
  if (!is.na(x$model)){
    cat(paste0("\tRequested model: ",x$model))
  }
  cat("\n")
  if (x$xreg.n > 0){
    if (x$pr.comp > 0){
      cat(paste0("Xreg: ",x$xreg.n," variables, reduced to ",x$pr.comp," principal components\n"))
    } else {
      cat(paste0("Xreg: ",x$xreg.n," variables\n"))
    }
  }
  cat("\n")
  
  tbl <- x$table
  if (all(tbl$Xreg == 0)){
    tbl$Xreg <- NULL
  }
  # Drop parameters that are not used at any aggregation level
  for (prmname in c("beta","gamma","phi")){
    if (all(is.na(tbl[,prmname]))){
      tbl[,prmname] <- NULL
    }
  }
  num <- unlist(lapply(tbl,is.double))
  tbl[num] <- lapply(tbl[num],function(v){round(v,digits)})
  print(tbl, row.names=FALSE, na.print="")
  cat("n: aggregated sample size; m: seasonal periods at aggregation level.\n")
  if (any(tbl$Model == "-")){
    cat("Levels marked '-' were not estimated due to insufficient sample.\n")
  }
  
  return(invisible(x))
  
}

#-------------------------------------------------
#' Produce estimated MAPA fit plot
#'
#' @description
#' The following function plots the fitted components at each aggregation
#' level of an estimated MAPA. The plot is the same as the one that can be produced optionally by mapaest. Aggregation levels for which seasonality cannot be estimated are in grey and levels that were not estimated due to insufficient sample are marked with '-'.
#'
#' @param x Fitted MAPA model (from mapaest).
#' @param xreg.plot Add infromation about xreg in the figure.
#' @param col Vector of 5 colours for: none, additive, additive damped, multiplicative, and multiplicative damped. If NULL the default colours are used.
#' @param ... Additional arguments passed to \code{image}, for example \code{main}, \code{xlab}, and \code{ylab}.
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @seealso \code{\link{mapaest}}, \code{\link{mapa}}.
#'
#' @keywords ~plotmapa ~plot
#' @examples
#' mapafit <- mapaest(admissions)
#' plot(mapafit)
#' plot(mapafit, main="Admissions")
#'
#' # With es and xreg, the number of xreg variables is shown at each level
#' set.seed(1)
#' x <- rnorm(length(admissions))
#' mapafit.x <- mapaest(admissions, type="es", xreg=x)
#' plot(mapafit.x)
#' plot(mapafit.x, xreg.plot=FALSE, col=heat.colors(5))
#' @export
plot.mapa.fit <- function(x,xreg.plot=c(TRUE,FALSE),col=NULL,...){
  xreg.plot <- xreg.plot[1]
  
  # Find locations in mapafit
  x <- mapafit.convert(x)
  
  # Get settings from mapafit
  ALs <- x$AL
  ets.type <- x$type
  ppy <- x$ppy
  k <- length(ALs)
  
  # Get permitted seasonal levels
  perm_seas <- (ppy %% ALs) == 0 & (ALs < ppy)
  
  # Fitted components at each AL, coded as numbers 
  # Levels that are not estimated are NA
  comps <- array(NA,c(k,3))
  comps.char <- array("",c(k,3))
  x.n <- vector("numeric",k)
  
  for (i in which(x$use)){
    
    # Get fitted components at each AL
    components <- mapa.components(x$models[[i]],ets.type)
    
    # Check for xreg
    if (ets.type == "es"){
      x.n[i] <- length(x$models[[i]]$initial$xreg)
    }
    
    # Error term
    if (components[1]=="A"){
      comps[i,1] <- 1
    } else {
      comps[i,1] <- 2
    }
    # Trend term
    comps[i,2] <- switch(components[2], "A"=1, "Ad"=1.5, "M"=2, "Md"=2.5, 0)
    # Season term
    comps[i,3] <- switch(components[3], "A"=1, "M"=2, 0)
    comps.char[i,] <- components
  }
  ttl <- paste(toupper(ets.type), "components")
  
  # Add row for xreg if needed
  if (sum(x.n)>0 & xreg.plot == TRUE){
    x.add <- 1
    comps <- cbind(comps,x.n>0)
    comps[!x$use,4] <- NA
  } else {
    x.add <- 0
  }
  nr <- 3 + x.add
  
  if (is.null(col)){
    col <- brewer.pal(8,"PuBu")[2:6] #rev(heat.colors(5)),
  }
  
  # User arguments override defaults
  args <- list(x=ALs, y=1:nr, z=matrix(comps[,nr:1],ncol=nr), axes=FALSE, col=col, 
               ylab="Components", xlab="Aggregation Level", main=ttl, 
               breaks=c(-1,0,1,1.5,2,2.5))
  args <- modifyList(args,list(...))
  do.call(image,args)
  axis(2, at=1:nr, labels=list("Error","Trend","Season","Xreg")[nr:1])
  axis(1, at=ALs)
  # Grey seasonality
  kk <- 1+x.add 
  for (AL in ALs[!perm_seas]){
    polygon(c((AL-0.5)*c(1,1),(AL+0.5)*c(1,1)),c(.5+kk,-.5+kk,-.5+kk,.5+kk),col="gray60",border=NA)
  }
  box()
  
  # Grid and labels
  for (AL in ALs){
    lines(c(AL-0.5,AL-0.5),c(0,nr+1),col="black")
  }
  for (i in 1:nr){
    lines(c(min(ALs)-0.5,max(ALs)+0.5),c(i-0.5,i-0.5),col="black")
  }
  for (i in 1:k){
    if (x$use[i]){
      text(ALs[i],(1+x.add):nr,rev(comps.char[i,]))
      if (x.add==1){
        text(ALs[i],1,x.n[i])
      }
    } else {
      text(ALs[i],1:nr,"-")
    }
  }
  
  return(invisible())

}
