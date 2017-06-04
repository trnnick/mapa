# These functions handle forecasting using MAPA
# Nikolaos Kourentzes 2016

statetranslate <- function(fit,AL,fh,q,ppyA,fittype,ets.type){
  # This function prepares ets states for MAPA combination
  # It extrapolates from last states the forecasts and translates to additive
  
  FCs_temp <- array(0, c(5, fh))
  fhA <- ceiling(fh/AL)   # forecast horizon for the aggregated level
  
  if (ets.type == "ets"){
    # Handle ets from forecast package
  
    param <- fit$par
    pnames <- names(param)
    if (sum(pnames=="phi")>0) {
      phi <- param[pnames=="phi"]
    } else {
      phi <- NULL
    }
    
    # Estimates for the Level Component
    FCs_temp[2, ] <- as.numeric(rep(rep(fit$states[q+1, 1], fhA), each=AL)[1:fh])
    
    # Estimates for the Trend Component
    if (fit$components[2]=="N"){ # no trend
      FCs_temp[3, ] <- 0
      b = 0 # indicates that there is no trend
    } else if (fit$components[2]=="A"){ # additive trend
      if (fit$components[4]=="FALSE"){ 
        FCs_temp[3, ] <- as.numeric(rep(fit$states[q+1, 2] * (1:fhA), each=AL))[1:fh]
      } else { # additive damped trend
        # We divide with phi because of an internal calculation for 
        # the damped trend in the ETS package
        FCs_temp[3, ] <- as.numeric(rep(cumsum((fit$states[q+1, 2]/phi)* phi^(1:fhA)), each=AL))[1:fh]
      }
      b <- 1 # indicates that there is trend
    } else {
      if (fit$components[4]=="FALSE"){ # multiplicative trend
        FCs_temp[3, ] <- as.numeric(rep((fit$states[q+1,2]^(1:fhA)-1), each=AL)[1:fh] * FCs_temp[2,])
      } else { # multiplicative damped trend
        # We divide with phi because of an internal calculation for the damped trend in the ETS package
        FCs_temp[3, ] <- as.numeric(rep((((fit$states[q+1,2] ^ (1/phi)) ^ cumsum(phi^(1:fhA)))-1), 
                                        each=AL)[1:fh] * FCs_temp[2, ])
      }
      b <- 1 # indicates that there is trend
    }
    
    # Estimates for the Seasonal Component 
    if (fit$components[3]=="N"){ # no seasonality
      FCs_temp[4, ] <- 0
    } else if (fit$components[3]=="A"){ # additive seasonality
      FCs_temp[4, ] <- as.numeric(rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA), each=AL))[1:fh]
    } else { # multiplicative seasonality
      FCs_temp[4, ] <- as.numeric((rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA),
                                       each=AL)[1:fh] - 1)) * (FCs_temp[2, ] + FCs_temp[3, ])
    }
    
  } else if (ets.type == "es"){
    # Handle es from smooth package
    
    init.obs <- max(1,length(fit$initialSeason))
    phi <- fit$phi
    
    # Get model components
    components <- fit$model
    components <- substring(components,5,nchar(components))
    components <- substring(components,1,nchar(components)-1)
    mn <- nchar(components)
    components <- substring(components, seq(1,mn,1), seq(1,mn,1))
    if (mn == 4){
      components <- c(components[1],paste0(components[2],components[3]),components[mn])
    }
    
    # Estimates for the Level Component - Careful es updates future level!!!
    FCs_temp[2, ] <- as.numeric(rep(rep(fit$states[q+init.obs, 1], fhA), each=AL)[1:fh])
    
    # Estimates for the Trend Component
    if (components[2]=="A"){ # additive trend
      FCs_temp[3, ] <- as.numeric(rep(fit$states[q+1+init.obs, 2] * (1:fhA), each=AL))[1:fh]
    } else if (components[2]=="Ad"){ # additive damped trend
      FCs_temp[3, ] <- as.numeric(rep(cumsum((fit$states[q+1+init.obs, 2]/phi)* phi^(1:fhA)), each=AL))[1:fh]
    } else if (components[2]=="M"){ # multiplicative trend
      FCs_temp[3, ] <- as.numeric(rep((fit$states[q+1+init.obs,2]^(1:fhA)-1), each=AL)[1:fh] * FCs_temp[2,])
    } else if (components[2]=="Md"){ # multiplicative damped trend
      FCs_temp[3, ] <- as.numeric(rep((((fit$states[q+1+init.obs,2] ^ (1/phi)) ^ cumsum(phi^(1:fhA)))-1), 
                                      each=AL)[1:fh] * FCs_temp[2, ])  
    } else {
      # No trend
      FCs_temp[3, ] <- 0
    }
    
    # Estimates for the Seasonal Component 
    # sc <- dim(fit$states)[2]
    sc <- which(colnames(fit$states) == "seasonality")
    if (components[3]=="N"){ # no seasonality
      FCs_temp[4, ] <- 0
    } else if (components[3]=="A"){ # additive seasonality
      FCs_temp[4, ] <- as.numeric(rep(rep(fit$states[(q+1+init.obs):(q+ppyA+init.obs),sc], fhA), each=AL))[1:fh]
    } else { # multiplicative seasonality
      FCs_temp[4, ] <- as.numeric((rep(rep(fit$states[(q+1+init.obs):(q+ppyA+init.obs),sc], fhA),
                                       each=AL)[1:fh] - 1)) * (FCs_temp[2, ] + FCs_temp[3, ])
    }
    
    # Estimate for the xreg 
    xreg <- fit$xreg
    if (!is.null(xreg)){
      if (!is.null(dim(xreg))){
        FCs_temp[5, ] <- as.numeric(rep( colSums(t(xreg[(q+1):(q+fhA),] * matrix(rep(fit$initialX,fhA),nrow=fhA,byrow=TRUE))) , each=AL))[1:fh]
      } else {
        FCs_temp[5, ] <- as.numeric(rep(xreg[(q+1):(q+fhA)] * fit$initialX, each=AL))[1:fh]
      }
    }
    
  }
  
  # fittype identifies if information is comming from ets or mapafit
  if (fittype==1){
    # Recreate ETS forecasts
    if (fh != 1) {
      FCs_temp[1, ] <- colSums(FCs_temp[2:5,])     
    } else {
      FCs_temp[1, ] <- sum(FCs_temp[2:5,])   
    }
  }
  
  # Return output
  return(FCs_temp)
  
}

#-------------------------------------------------
mapacalc <- function(y, mapafit, fh=0, comb=c("mean", "median", "w.mean", "w.median", "wght"),
                     outplot=c(0,1,2), hybrid=c(TRUE,FALSE), xreg=NULL){
  # Calculation of MAPA forecasts
  # 
  # Inputs:
  #   y           = In sample observations of a time series (vector)
  #   mapafit     = Fitted MAPA model (from mapaest)
  #   fh          = Forecast horizon. Default = ppy
  #   comb        = Combination operator. This can be: "mean"; "median"; "wght" - where each 
  #                 aggregation level is weighted inversly to aggregation; "w.mean" - level 
  #                 and trend components are averaged, but seasonal and xreg follow the wght 
  #                 combination; "w.median" - as w.mean, but with median. It is suggested that 
  #                 for data with high sampling frequency to use one of the "w.mean" and "w.median".
  #   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
  #                 2 = time series, forecasts and components. Default is 1. 
  #   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
  #                 If minimumAL > 1 then the minimumAL ETS forecasts are used.
  #   xreg        = Vector or matrix of exogenous variables to be included in the MAPA. 
  #                 If matrix then rows are observations and columns are variables. 
  #                 Must be at least as long as in-sample plus fh. Additional observations are unused.
  #
  # Output:
  #   forecasts   = Vector with forecasts
  #   components  = array with MAPA components
  
  comb <- match.arg(comb,c("mean","median","w.mean","w.median","wght"))
  outplot <- outplot[1]
  hybrid <- hybrid[1]
  
  # Find locations in mapafit
  idx.use <- which(colnames(mapafit)=="use")
  idx.AL <- which(colnames(mapafit)=="AL")
  idx.ppy <- which(colnames(mapafit)=="original.ppy")
  idx.comp <- which(colnames(mapafit)=="pr.comp")
  ets.type <- mapafit[1,which(colnames(mapafit)=="etstype")]
  
  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,idx.use]==TRUE, idx.AL])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(mapafit[1,idx.ppy])
  
  # Set default foreast horizon
  if (fh == 0){
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
    if (!is.null(dim(xreg))){
      p <- dim(xreg)[2]
    } else {
      p <- 1
    }
    xreg <- matrix(xreg,ncol=p)
    # Check that it contains at least n observatoins and trim
    xn <- dim(xreg)[1]
    if (xn < (fh + observations)){
      stop("Number of observations of xreg must be equal of exceed the length of y + fh.")
    } 
  }
  
  # The forecasted components are saved here
  FCs <- array(0, c(maximumAL-minimumAL+1, 5, fh),dimnames=list(paste("AL",minimumAL:maximumAL,sep=""),
                                                                c("ETS","Level","Trend","Season","Xreg"),paste("t+",1:fh,sep=""))) 
  
  # MAPA forecast
  ALvec <- minimumAL:maximumAL
  
  for (ALi in 1:(maximumAL-minimumAL+1)){
    
    AL <- ALvec[ALi]
    fhA <- ceiling(fh/AL)   # This is used for es and xreg preprocessing
    
    q <- observations %/% AL # observation in the aggregated level
    ppyA <- ppy %/% AL       # periods per year for the aggregated level
    if (ppy %% AL != 0){
      ppyA <- 1
    }
    
    # Aggregation
    yA <- colMeans(matrix(tail(y,q*AL),nrow=AL))
    ats <- ts(yA, frequency = ppyA) 
    
    # Aggregate xreg 
    # This is aligned with the end point of the aggregation of y
    if (!is.null(xreg)){
      r <- observations - q*AL
      p <- dim(xreg)[2]
      xn <- dim(xreg)[1]
      xregA <- array(NA, c(q+fhA,p))
      for (k in 1:p){
        temp <- xreg[(r+1):min(c((q+fhA)*AL,xn)),k]
        m <- ceiling(length(temp)/AL)*AL - length(temp)
        temp <- c(temp,rep(NA,m))
        xregA[,k] <- colMeans(matrix(temp,nrow=AL),na.rm=TRUE)
      }

      # Check that mapafit expected xreg matches xreg input
      # Either pr.comp = 0 and check number of x's or pr.comp > 0 and x's should be >= pr.comp
      idx.x <- which(colnames(mapafit)=="initialX")
      if (((mapafit[[ALi,idx.comp]]$pr.comp > 0) & (length(mapafit[[ALi,idx.comp]]$mean) != p)) |
          ((mapafit[[ALi,idx.comp]]$pr.comp == 0) & (length(mapafit[[ALi,idx.x]]) != p))){
        stop("Number of xreg input variables does not match mapafit specification.")
      }
      # Preprocess if needed
      pr.out <- mapaprcomp(xregA,mapafit[[ALi,idx.comp]])
      xregA <- pr.out$x.out
      
    } else {
      xregA <- NULL
    }
    
    # Predict
    if (ets.type == "ets"){
      
      AL.fit <- structure(mapafit[ALi,1:18],class="ets")
      ats.fit <- ets(ats, AL.fit, use.initial.values=TRUE)

    } else if (ets.type=="es"){
      
      AL.fit <- structure(mapafit[ALi,1:32], class = "smooth")
      ats.fit <- suppressWarnings(es(ats, model=AL.fit, h=fhA, silent="all",xreg=xregA))
      
    }
    
    # Transalte ets states for MAPA
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,1,ets.type)
    
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
  combres <- mapacomb(minimumAL,maximumAL,ppy,FCs,comb)
  forecasts <- combres$forecasts
  perm_levels <- combres$perm_levels
  perm_seas <- combres$perm_seas
  
  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }  

  # Plot output
  mapaplot(outplot,FCs,minimumAL,maximumAL,perm_levels,perm_seas,observations,y,forecasts,fh,comb)
  
  # Construct output
  return(structure(list("forecast"=forecasts,"components"=FCs),class="mapa.calc"))
  
}

#-------------------------------------------------
summary.mapa.calc <- function(object,...){
  print(object)
}

#-------------------------------------------------
print.mapa.calc <- function(x,...){
  print(x$forecast)
}

#-------------------------------------------------
mapafor <- function(y, mapafit, fh=-1, ifh=1, comb=c("mean","median","w.mean","w.median","wght"), 
                    outplot=c(1,0), hybrid=c(TRUE,FALSE), conf.lvl=NULL, xreg=NULL) {
  # MAPA in- and out-of-sample forecast
  # 
  # Inputs:
  #   y           = In sample observations of a time series (vector)
  #   mapafit     = Fitted MAPA model (from mapaest)
  #   fh          = Forecast horizon. Default = ppy
  #   ifh         = In-sample forecast horizon. Default = 1
  #   comb        = Combination operator. This can be: "mean"; "median"; "wght" - where each 
  #                 aggregation level is weighted inversly to aggregation; "w.mean" - level 
  #                 and trend components are averaged, but seasonal and xreg follow the wght 
  #                 combination; "w.median" - as w.mean, but with median. It is suggested that 
  #                 for data with high sampling frequency to use one of the "w.mean" and "w.median".
  #   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1. 
  #   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
  #                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
  #   conf.lvl    = Vector of confidence level for prediction intervals. Values must be (0,1). 
  #                 If conf.lvl == NULL then no intervals are calculated. For example to get 
  #                 the intervals for 80% and 95% use conf.lvl=c(0.8,0.95). 
  #   xreg        = Vector or matrix of exogenous variables to be included in the MAPA. 
  #                 If matrix then rows are observations and columns are variables. 
  #                 Must be at least as long as in-sample plus fh. Additional observations are unused.
  #
  # Output:
  #   out$infor   = In-sample forecasts.
  #   out$outfor  = Out-of-sample forecasts.
  #   out$PI      = Prediction intervals for given confidence levels.
  #   out$MSE     = In-sample MSE error.
  #   out$MAE     = In-sample MAE error.
  
  comb <- match.arg(comb,c("mean","median","w.mean","w.median","wght"))
  outplot <- outplot[1]
  hybrid <- hybrid[1]
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Find locations in mapafit
  idx.use <- which(colnames(mapafit)=="use")
  idx.AL <- which(colnames(mapafit)=="AL")
  idx.ppy <- which(colnames(mapafit)=="original.ppy")
  ets.type <- mapafit[1,which(colnames(mapafit)=="etstype")]
  
  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,idx.use]==TRUE, idx.AL])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  
  ppy <- as.numeric(mapafit[[1,idx.ppy]])
  
  if (fh == -1){
    fh <- ppy
  }
  
  # Override ifh.c to be >= fh if conf.lvl is not NULL
  # Output will crop values to ifh
  ifh.c <- ifh
  if (!is.null(conf.lvl)){
    if (ifh.c < fh){
      ifh.c <- fh
    } 
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
    for (i in i.start:(observations-1)){
      inobs <- as.matrix(y[1:i])
      infor[, i+1] <- mapacalc(inobs, mapafit, fh=ifh.c, comb, outplot=0, hybrid, xreg=xreg)$forecast
      # Crop out-of-sample predictions
      if ((i+ifh.c)>observations){
        k <- (i+ifh.c) - observations
        infor[(ifh.c-k+1):ifh.c, i+1] <- rep(NA,k)
      }
    }    
  } else {
    infor <- NULL
    ifh.c <- 0        # Override in-sample output if no forecasts are calculated
    conf.lvl <- NULL  # Confidence intervals cannot be calculated
    
  }
  
  # Out-of-sample MAPA
  if (fh>0){
    outfor <- mapacalc(y, mapafit, fh, comb, outplot=0, hybrid,xreg=xreg)$forecast
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
  
  # Produce plot
  if (outplot==1){
    layout(matrix(1, 1, 1, byrow = TRUE))
    # Find min max
    if (is.null(outfor)){
      ymax <- max(y)
      ymin <- min(y)
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)      
    } else {
      if (!is.null(conf.lvl)){
        ymax <- max(c(max(outfor),max(y),max(PI)))
        ymin <- min(c(min(outfor),min(y),min(PI)))
      } else {
        ymax <- max(c(max(outfor),max(y)))
        ymin <- min(c(min(outfor),min(y)))
      }
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)
    }
    cmp <- brewer.pal(3,"Set1")
    plot(1:observations,y,type="l",col=cmp[2], xlab="", ylab="", main="Forecast", 
         xlim <- c(1, observations+fh), ylim=c(ymin,ymax))
    # In-sample
    if (ifh.c>0){
      if (ifh==1){
        lines(infor[1,],col=cmp[1])
      } else {
        # clrs = rainbow(observations-ppy)
        for (i in (ppy):(observations-1)){
          lines((i):(i+ifh-1),infor[,i],col=cmp[1])
        }
      }
    }
    # Prediction intervals
    if (!is.null(conf.lvl)){
      cmp2 <- rgb(1,0.75,0.75,1)
      for (i in 1:PIn){
        rr <- 1-(1-(i-1)*(0.5/PIn))*(1-c(col2rgb(cmp2)/255))
        polygon(c(observations+(1:fh),observations+(fh:1)),c(PI[i,],PI[PIn*2+1-i,fh:1]),
                col=rgb(rr[1],rr[2],rr[3]), border=NA)
      }
    }

    # Out-of-sample
    if (ifh == 0 || ifh.c == 0){
      lines(observations:(fh+observations),c(y[observations],outfor),col=cmp[1])
    } else if (ifh == 1){
      lines(observations:(fh+observations),c(infor[1,observations],outfor),col=cmp[1])
    } else {
      lines((observations+1):(fh+observations),outfor,col=cmp[1])
    }
  }
  
  # Construct output
  return(structure(list("infor"=infor,"outfor"=outfor,"PI"=PI,"MSE"=MSE,"MAE"=MAE),class="mapa.for"))
  
}

#-------------------------------------------------
summary.mapa.for <- function(object,...){
  print(object)
}

#-------------------------------------------------
print.mapa.for <- function(x,...){
  
  if (is.null(x$MSE)){
    print(x$outfor)
  } else {
    cat(paste0("MAPA fit MSE: ",round(x$MSE,2), ", MAE: ", round(x$MAE,2),"\n"))
    cat("Out-of-samplpe forecasts:\n")
    if (is.null(x$PI)){
      print(x$outfor)
    } else {
      sz <- dim(x$PI)[1]
      temp <- t(rbind(x$PI[(1+sz/2):sz,],x$outfor,x$PI[1:(sz/2),]))
      colnames(temp)[sz/2 + 1] <- "Forecast"
      print(temp)
    }
  }
  
}

