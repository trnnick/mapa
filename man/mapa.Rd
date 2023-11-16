\name{mapa}
\alias{mapa}

\title{Multiple Aggregation Prediction Algorithm (Wrapper)}

\description{Wrapper to estimate and produce MAPA in- and out-of-sample forecasts. 
Uses mapaest and mapafor.}

\usage{
mapa(y, ppy, fh=ppy, ifh=1, minimumAL=1, maximumAL=ppy, 
     comb=c("w.mean","w.median","mean","median","wght"), paral=c(0,1,2), 
     display=c(0,1), outplot=c(0,1), hybrid=c(TRUE,FALSE), model="ZZZ", 
     type=c("ets","es"), conf.lvl=NULL, xreg=NULL, pr.comp=0, ...)
}

\arguments{
  \item{y}{
    In sample observations of a time series (vector). If y == "paper" then it prints paper reference.
}
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{ifh}{
    Lower aggregation level to use. Default = 1.
}
  \item{minimumAL}{
    Lowest aggregation level to use. Default = 1.
}
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy, maximumAL>1.
}
  \item{comb}{
    Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
}
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes and initialise cluster. Default is 0.
}
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 1.
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{model}{
    Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M" or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X" or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X" or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. 
}
  \item{type}{
    What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
}
  \item{conf.lvl}{
    Vector of confidence level for prediction intervals. Values must be (0,1). If conf.lvl == NULL then no intervals are calculated. For example to get the intervals for 80\% and 95\% use conf.lvl=c(0.8,0.95).
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
\item{infor}{In-sample forecasts.}
\item{outfor}{Out-of-sample forecasts.}
\item{PI}{Prediction intervals for given confidence levels.}
\item{MSE}{In-sample MSE error.}
\item{MAE}{In-sample MAE error.}
}
\details{This function calls \emph{mapaest} and \emph{mapafor} internally.}
\note{The calculation of the prediction intervals is based on the empirical multiple step ahead MSE. To speed up calculations set \code{conf.lvl=NULL}. If very long forecast horizons are requested then once no more t+h MSE can be calculated the following approximation is used: \eqn{\sqrt{MSE_{t+1}}*\sqrt{h}}{sqrt(MSE_[t+1])*sqrt(h)} for the error.}
\references{
    \itemize{
        \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time 
series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
        \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
        \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
}

}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}; Fotios Petropoulos.
}
\seealso{
\code{\link{mapaest}}, \code{\link{mapafor}}, \code{\link{mapasimple}}.
}
\examples{
out <- mapa(admissions)
}
\keyword{ ~mapaest }
\keyword{ ~mapafor }
\keyword{ ~mapasimple }