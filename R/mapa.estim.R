# These functions are handling the estimation of MAPA
# Nikolaos Kourentzes 2016 

mapaest <- function(y, ppy=NULL, minimumAL=1, maximumAL=ppy, paral=c(0,1,2), display=c(0,1), 
                    outplot=c(0,1), model="ZZZ", type=c("ets","es"), ...) {
  # Estimate MAPA for a time series  
  #  
  # Inputs:
  #   y           = In sample observations of a time series (vector)
  #   ppy         = Periods in a season of the time series at the sampled frequency.
  #                 If y is a ts object then this is taken from its frequency,
  #                 unless overriden. 
  #   minimumAL   = Lowest aggregation level to use. Default = 1, maximumAL>1
  #   maximumAL   = Highest aggregation level to use. Default = ppy
  #   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
  #                 2 = yes and initialise cluster. Default is 0.
  #   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
  #   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.  
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
  #   ...         = Pass additional arguments to es and ets.
  #
  # Output:
  #   mapafit     = Estimated MAPA model structure
  
  # Default type
  paral <- paral[1]
  display <- display[1]
  outplot <- outplot[1]
  type <- type[1]
  
  # Get ppy and maximumAL
  if (is.null(ppy)){
    if (class(y)=="ts"){
      ppy <- frequency(y)
      if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
      stop(paste("Input ppy is not given and y input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }  
  
  # Suggest to use es instead of ets when ppy > 24
  if (ppy > 24 & type == "ets"){
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
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
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
  if (paral != 0){  # Parallel run
    mapafit <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapaest.loop, 
                              y=y, minimumAL=minimumAL, maximumAL=maximumAL, observations=observations, ppy=ppy,
                              display=display,model=model,type=type,...)  
  } else {          # Serial run
    mapafit <- vector("list", (maximumAL-minimumAL+1))
    for (i in 1:(maximumAL-minimumAL+1)){
      mapafit[[i]] <- mapaest.loop(i, y, minimumAL, maximumAL, observations, 
                                   ppy, display,model=model,type=type,...)
    }
  }
  
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }
  
  # Process output
  mapafit <- do.call(rbind, mapafit) # Re-arrange output for clusterApplyLB function
  rownames(mapafit) <- paste("AL",minimumAL:maximumAL,sep="")
  mapafit <- structure(mapafit, class = "mapa.fit")
  
  # Plot model selection summary
  if (outplot == 1){
    plot(mapafit)
  }
    
  # Return output
  return(mapafit)

}

#-------------------------------------------------
mapaest.loop <- function(ALi, y, minimumAL, maximumAL, observations, 
                         ppy, display, model, type, ...){ 
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
  
  q <- observations %/% AL # observation in the aggregated level
  # r <- observations %% AL  # observation to discard from the beginning of the series
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
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
  
  if (q >= 4){ # Fit only when there is enough sample
    
    # Aggregation
    yA <- colMeans(matrix(tail(y,q*AL),nrow=AL))
    ats <- ts(yA, frequency = ppyA) # create the time series
    # Old code
    # yA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    # for (j in 1:q){                 # calculate the aggregate values
    #   yA[j] <- mean(y[(r+1+(j-1)*AL):(r+j*AL)])
    # }
    
    # Check if seasonality exists and select appropriate model
    if ((tail(model.set,1) == "A" | tail(model.set,1) == "M") && ppyA == 1){
      mapa.model <- paste0(model.set[1],model.set[2],"N")
    } else if ((ppyA == 2 && q <= 3*ppyA) & type == "ets"){   
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
      # structure, correct for that
      if (length(names(fit)) != 18){
        fit.temp = fit
        fit <- list("loglik"=NULL,"aic"=NULL,"bic"=NULL,"aicc"=NULL,"mse"=NULL,
                    "amse"=NULL,"fit"=NULL,"residuals"=NULL,"fitted"=NULL,
                    "states"=NULL,"par"=NULL,"m"=NULL,"method"=NULL,
                    "components"=NULL,"call"=NULL,"initstate"=NULL,
                    "sigma2"=NULL,"x"=NULL)
        fit.temp.names <- names(fit.temp)
        fit.names <- names(fit)
        for (fi in 1:length(fit.names)){
          if (sum(fit.temp.names == fit.names[fi])>0){
            eval(parse(text=paste0("fit$",fit.names[fi]," <- fit.temp$",fit.names[fi])))
          }
        }
      }
      
    } else if (type == "es"){
      
      # Turn off warnings for es - this is done when the model reduces pool 
      # due to sample size.
      fit <- suppressWarnings(es(ats,model=mapa.model,silent="all", ...))
            
    }
    
    # Common part
    fit$use <- TRUE
    
  } else {
    # No fit is taking place
    
    if (type == "ets"){
      fit <- list("loglik"=NULL,"aic"=NULL,"bic"=NULL,"aicc"=NULL,"mse"=NULL,
                  "amse"=NULL,"fit"=NULL,"residuals"=NULL,"fitted"=NULL,
                  "states"=NULL,"par"=NULL,"m"=NULL,"method"=NULL,
                  "components"=NULL,"call"=NULL,"initstate"=NULL,
                  "sigma2"=NULL,"x"=NULL,"use"=FALSE)
    } else if (type == "es"){
      fit <- list("model"=NULL,"timeElapsed"=NULL,"states"=NULL,
                  "persistence"=NULL,"phi"=NULL,"initialType"=NULL,
                  "initial"=NULL,"initialSeason"=NULL,"nParam"=NULL,
                  "fitted"=NULL,"forecast"=NULL,"lower"=NULL,
                  "upper"=NULL,"residuals"=NULL,"errors"=NULL,
                  "s2"=NULL,"intervals"=NULL,"level"=NULL,
                  "actuals"=NULL,"holdout"=NULL,"iprob"=NULL,
                  "intermittent"=NULL,"xreg"=NULL,"updateX"=NULL,
                  "initialX"=NULL,"persistenceX"=NULL,"transitionX"=NULL,
                  "ICs"=NULL,"cf"=NULL,"cfType"=NULL,"FI"=NULL,
                  "accuracy"=NULL,"use"=FALSE)
    }
    
  }
  
  fit$AL <- AL
  fit$original.ppy <- ppy
  if (type == "ets"){
    fit$etstype <- "ets"
  } else if (type == "es"){
    fit$etstype <- "es"
  }
  
  # Update console display
  if (display==1){
    nc = nchar(txtc)
    cat(rep("\r",nc))
    cat(rep(" ",nc))
    cat(rep("\r",nc))
  }
  
  # Return loop result
  return(rbind(fit))
  
}

#-------------------------------------------------
plotmapa <- function(mapafit){
  
  warning("This function is depreciated. Simply use plot().")
  plot.mapa.fit(mapafit)
  
}

#-------------------------------------------------
summary.mapa.fit <- function(object,...){
  print(object)
}

#-------------------------------------------------
print.mapa.fit <- function(x,...){
  
  # Get basic settings
  idx.use <- which(colnames(x)=="use")
  idx.AL <- which(colnames(x)=="AL")
  idx.ppy <- which(colnames(x)=="original.ppy")
  ets.type <- x[1,which(colnames(x)=="etstype")]
  
  # Get settings from mapafit
  ALs <- as.numeric(x[x[,idx.use]==TRUE, idx.AL])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(x[1,idx.ppy])
  
  # Plot model selection summary
  ALplot <- 1:(maximumAL-minimumAL+1)
  ALplot <- ALplot[unlist(x[,idx.use])==TRUE]
  
  cat(paste0("MAPA fitted using ",ets.type, "\tOriginal frequency: ",ppy,"\n"))
  if (ets.type == "ets"){
    mdls <- x[unlist(x[,idx.use]),13]
    mdls <- lapply(mdls,function(x){gsub(',','',x)})
  } else if (ets.type == "es"){
    mdls <- x[unlist(x[,idx.use]),1]
  }
  for (i in 1:(maximumAL-minimumAL+1)){
    cat(paste0("Aggregation level: ",ALplot[i],"\t","Method: ",mdls[[i]],"\n"))
  }
    
}

#-------------------------------------------------
plot.mapa.fit <- function(x,...){
  # Produce estimated MAPA fit plot
  # 
  # Inputs:
  #   mapafit     = Fitted MAPA model (from mapaest)
  
  # Find locations in mapafit
  idx.use <- which(colnames(x)=="use")
  idx.AL <- which(colnames(x)=="AL")
  idx.ppy <- which(colnames(x)=="original.ppy")
  ets.type <- x[1,which(colnames(x)=="etstype")]
  
  # Get settings from mapafit
  ALs <- as.numeric(x[x[,idx.use]==TRUE, idx.AL])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(x[1,idx.ppy])
  
  # Plot model selection summary
  ALplot <- 1:(maximumAL-minimumAL+1)
  ALplot <- ALplot[unlist(x[,idx.use])==TRUE]
  
  layout(matrix(1, 1, 1, byrow = TRUE))
  comps <- array(0,c(max(ALplot),4))
  comps.char <- array(0,c(max(ALplot),3))
  for (AL in 1:max(ALplot)){
    
    # Get fitted components at each AL
    if (ets.type == "ets"){
      components <- x[[AL, 14]]
      if (components[[4]] == "TRUE"){
        components[[4]] <- "d"
      } else {
        components[[4]] <- ""
      }
      components <- components[c(1,2,4,3)]
      components <- c(components[1],paste0(components[2],components[3]),tail(components,1))
      ttl <- "ETS components"
    } else if (ets.type == "es"){
      # Get ES components
      components <- x[[AL, 1]]
      components <- substring(components,5,nchar(components))
      components <- substring(components,1,nchar(components)-1)
      mn <- nchar(components)
      components <- substring(components, seq(1,mn,1), seq(1,mn,1))
      if (mn == 4){
        components <- c(components[1],paste0(components[2],components[3]),components[mn])
      }
      ttl <- "ES components"
    }
    
    # Error term
    if (components[1]=="A"){
      comps[AL,1] <- 1
    } else {
      comps[AL,1] <- 2
    }
    comps.char[AL,1] <- components[1]
    # Trend term
    if (components[2]=="A"){
      comps[AL,2] <- 1
    } else if (components[2]=="Ad"){
      comps[AL,2] <- 1.5
    } else if (components[2]=="M"){
      comps[AL,2] <- 2
    } else if (components[2]=="Md"){
      comps[AL,2] <- 2.5
    } else {
      comps[AL,2] <- 0
    }
    comps.char[AL,2] <- components[2]
    # Season term
    if (components[3]=="A"){
      comps[AL,3] <- 1
    } else {if (components[3]=="M"){
      comps[AL,3] <- 2
    } else
      comps[AL,3] <- 0
    }
    comps.char[AL,3] <- components[3]
    comps[AL,4] <- x[[AL,idx.AL]]
  }
  
  image(min(comps[,4]):max(comps[,4]), 1:3, matrix(comps[,1:3],ncol=3), axes=FALSE, col=brewer.pal(8,"PuBu")[2:6], #rev(heat.colors(5)),
        ylab="Components", xlab="Aggregation Level", main=ttl,breaks=c(-1,0,1,1.5,2,2.5))
  axis(2, at=1:3, labels=list("Error","Trend","Season"))
  axis(1, at=min(comps[,4]):max(comps[,4]))
  box()
  
  for (i in 1:4){
    for (AL in 1:max(ALplot)){
      if (i==1){
        lines(c(AL-0.5+minimumAL-1,AL-0.5+minimumAL-1),c(0,4),col="black")
      }
      if (i<4 & AL<=max(comps[,4])){
        text(AL+minimumAL-1,i,comps.char[AL,i])
      }
    }
    lines(c(min(comps[,4])-0.5,max(comps[,4])+0.5),c(i-0.5,i-0.5),col="black")
  }

}
