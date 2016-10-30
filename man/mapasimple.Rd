\name{mapasimple}
\alias{mapasimple}
\title{MAPA Single Origin Forecast Wrapper}
\description{
The following function estimates MAPA and produces out-of-sample forecasts. 
It is not suggested, unless only simple single origin forecasts are required.
}
\usage{
mapasimple(y, ppy=NULL, fh=ppy, minimumAL=1, maximumAL=ppy, 
           comb=c("mean","median","wght"), paral=c(0,1,2), display=c(0,1), 
           outplot=c(1,0), hybrid=c(TRUE,FALSE), model="ZZZ", 
           type=c("ets","es"), xreg=NULL, pr.comp=0, ...)
}
\arguments{
  \item{y}{
    In sample observations of a time series (vector). If y == "paper" then 
	it prints paper reference.
}
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{minimumAL}{
    Lowest aggregation level to use. Default = 1, maximumAL>1.
}
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy.
}
  \item{comb}{
    Combination operator. One of "mean" or "median". Default is "mean". "wght" is experimental, use at your own risk!
}
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes 
	and initialise cluster. Default is 0.
}
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, forecasts and components. For the components the spectral colouring scheme is used. Dark red is aggregation level 1. Default is 1. 
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{model}{
    Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M" or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X" or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X" or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. For aggregation levels that seasonality becomes 1 then a non-seasonal model is estimated.
}
  \item{type}{
    What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
  }
  \item{xreg}{
    Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample. Additional observations are unused. Note that including xreg will force type="es". 
  }  
  \item{pr.comp}{
    MAPAx can use principal component analysis to preprocess xreg. When comp is -1 then the number of retained components is chosen automatically. When comp=0 then no pre-processing is performed and the original xreg is used. Any other value represents the number of principal components retained. 
  }
  \item{...}{
    Pass additional arguments to es or ets.
  }
}
\value{
\item{forecasts}{Vector with forecasts.}
\item{components}{Array with MAPA components, if output="all".}
}
\references{
Kourentzes N., Petropoulos F., Trapero J.R., 2014. Improving forecasting by estimating time 
series structural components across multiple frequencies, International Journal of Forecasting,
 30(2), 291-302.
}
\author{
Nikolaos Kourentzes and Fotios Petropoulos
}
\seealso{
\code{\link{mapa}}.
}
\examples{
mapasimple(admissions)
}
\keyword{ ~mapa }
