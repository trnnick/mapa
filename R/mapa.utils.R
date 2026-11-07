# Internal helper functions shared across MAPA estimation, forecasting, and plotting
# Nikolaos Kourentzes 2026

#-------------------------------------------------
mapa.paper <- function(){
# Print paper reference. Used when y is not numeric. 
  
  writeLines("Paper reference: ")
  writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
  writeLines("Improving forecasting by estimating time series structural components")
  writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
                   " 30 (2), 291-302.",sep=""))
  return(invisible())
  
}

#-------------------------------------------------
mapa.choice <- function(x, choices, name){
# Check that a numeric/logical option is one of the permitted values 
# and return the first element, which is the default.
  
  x <- x[1]
  if (is.na(x) || !(x %in% choices)){
    stop(paste0("Input ",name," must be one of: ",paste(choices,collapse=", "),"."))
  }
  return(x)
  
}

#-------------------------------------------------
mapa.alinfo <- function(AL, observations, ppy){
# Sample size and seasonal periods at an aggregation level
#
# Inputs:
#   AL            = Aggregation level.
#   observations  = Number of observations at the original sampling frequency.
#   ppy           = Periods in a season at the original sampling frequency.
#
# Output:
#   q             = Number of observations at the aggregation level.
#   ppyA          = Periods in a season at the aggregation level. If AL is not a 
#                   divisor of ppy then this is set to 1.
  
  q <- observations %/% AL # observation in the aggregated level
  # r <- observations %% AL  # observation to discard from the beginning of the series
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
  return(list("q"=q, "ppyA"=ppyA))
  
}

#-------------------------------------------------
mapa.aggr <- function(x, AL, q){
# Non-overlapping temporal aggregation, aligned with the end of the sample. 
# Any observations that do not fit a complete period are discarded from the 
# beginning of the series. 
#
# Inputs:
#   x             = Vector, or matrix with variables in columns.
#   AL            = Aggregation level.
#   q             = Number of aggregated observations to produce.
#
# Output:
#   xA            = Aggregated vector or matrix with q rows.
  
  if (is.null(dim(x))){
    xA <- colMeans(matrix(tail(x,q*AL),nrow=AL))
    # Old code
    # yA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    # for (j in 1:q){                 # calculate the aggregate values
    #   yA[j] <- mean(y[(r+1+(j-1)*AL):(r+j*AL)])
    # }
  } else {
    xA <- apply(x,2,function(xk){colMeans(matrix(tail(xk,q*AL),nrow=AL))})
    xA <- matrix(xA,nrow=q)
  }
  
  return(xA)
  
}

#-------------------------------------------------
mapa.components <- function(fit, ets.type){
# Get the components of a fitted exponential smoothing model 
# as c(error, trend, season), for example c("M","Ad","N").
  
  if (ets.type == "ets"){
    # Handle ets from forecast package
    components <- fit$components
    if (length(components) >= 4 && components[4] == "TRUE"){
      components <- c(components[1],paste0(components[2],"d"),components[3])
    } else {
      components <- components[1:3]
    }
  } else if (ets.type == "es"){
    # Handle es from smooth package
    # Get ES components
    components <- modelType(fit)
    mn <- nchar(components)
    components <- substring(components, seq(1,mn,1), seq(1,mn,1))
    if (mn == 4){
      components <- c(components[1],paste0(components[2],components[3]),components[mn])
    }
  }
  
  return(components)
  
}

#-------------------------------------------------
mapa.xreg <- function(xreg){
# Force xreg to be a matrix, each column a variable
  
  if (is.null(xreg)){
    return(NULL)
  }
  xreg <- as.matrix(xreg)
  if (!is.numeric(xreg)){
    stop("xreg must be numeric.")
  }
  
  return(xreg)
  
}

#-------------------------------------------------
mapa.cluster <- function(paral, cl=NULL){
# Prepare cluster for parallel processing
#
# Inputs:
#   paral         = 1 = use an existing cluster; 2 = initialise a new cluster.
#   cl            = An existing cluster, created with parallel::makeCluster.
#
# Output:
#   cl            = Cluster to use.
#   created       = TRUE if the cluster was initialised here and must be stopped.
  
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", length(cl), 'cores'))
    created <- TRUE
  } else {
    # For compatibility with earlier versions look for a cluster named cl, 
    # if one is not provided.
    if (is.null(cl)){
      cl <- get0("cl", envir=globalenv())
    }
    if (!inherits(cl,"cluster")){
      stop("paral=1 requires an initialised cluster. Provide it using the input cl, or use paral=2.")
    }
    created <- FALSE
  }
  
  return(list("cl"=cl, "created"=created))
  
}

#-------------------------------------------------
mapa.ylim <- function(...){
# Find min max for plot limits, with a 10% margin on either side
  
  yminmax <- range(unlist(list(...)), na.rm=TRUE, finite=TRUE)
  yminmax <- yminmax + c(-1,1)*0.1*diff(yminmax)
  
  return(yminmax)
  
}

#-------------------------------------------------
mapa.time <- function(y, fh=0){
# Time index of the in-sample observations (tin) and the next fh periods (tout).
# If y is a ts object its time attributes are used, otherwise 1, 2, ...
  
  n <- length(y)
  if (isa(y,"ts")){
    tin <- as.numeric(time(y))
    tout <- tin[n] + seq_len(fh)/frequency(y)
  } else {
    tin <- 1:n
    tout <- n + seq_len(fh)
  }
  
  return(list("tin"=tin, "tout"=tout))
  
}

#-------------------------------------------------
mapa.tsout <- function(x, y){
# If y is a ts object, return forecasts x as ts that start after the end of y
  
  if (isa(y,"ts") && !is.null(x)){
    x <- ts(as.numeric(x), start=tsp(y)[2]+1/frequency(y), frequency=frequency(y))
  }
  
  return(x)
  
}

#-------------------------------------------------
mapafit.convert <- function(mapafit){
# Convert mapafit objects created with MAPA <= 2.0.6 to the current structure. 
# Earlier versions stored the fitted models as a list-matrix, with one row 
# per aggregation level.
  
  if (!is.matrix(mapafit)){
    # Current structure, nothing to do
    return(mapafit)
  }
  
  cn <- colnames(mapafit)
  ets.type <- mapafit[[1,which(cn=="etstype")]]
  ppy <- as.numeric(mapafit[[1,which(cn=="original.ppy")]])
  pr.comp <- mapafit[[1,which(cn=="pr.comp")]]
  AL <- as.numeric(unlist(mapafit[,which(cn=="AL")]))
  use <- as.logical(unlist(mapafit[,which(cn=="use")]))
  k <- length(AL)
  
  models <- vector("list",k)
  q <- ppyA <- rep(NA,k)
  for (i in 1:k){
    ppyA[i] <- mapa.alinfo(AL[i],1,ppy)$ppyA
    if (use[i]){
      if (ets.type == "ets"){
        models[[i]] <- structure(mapafit[i,1:18],class="ets")
        q[i] <- length(models[[i]]$x)
      } else {
        smoothLength <- mapafit[[i,which(cn=="smoothLength")]]
        models[[i]] <- structure(mapafit[i,1:smoothLength],class=c("adam","smooth"))
        q[i] <- length(models[[i]]$fitted)
      }
    }
  }
  names(models) <- paste0("AL",AL)
  
  # Number of xreg variables
  if (pr.comp$pr.comp != 0){
    xreg.n <- length(pr.comp$mean)
  } else if (ets.type == "es"){
    xreg.n <- length(models[[which(use)[1]]]$initial$xreg)
  } else {
    xreg.n <- 0
  }
  
  mapafit <- structure(list("models"=models, "AL"=AL, "use"=use, "q"=q, "ppyA"=ppyA, 
                            "ppy"=ppy, "type"=ets.type, "model"=NA, "observations"=NA, 
                            "xreg.n"=xreg.n, "pr.comp"=pr.comp), class="mapa.fit")
  
  return(mapafit)
  
}
