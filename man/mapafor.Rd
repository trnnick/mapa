\name{mapafor}
\alias{mapafor}
\title{MAPA Forecast}
\description{
The following function produces in- and out-of-sample MAPA forecasts, for multiple steps ahead.
 This is the recommended function to use in forecasting with MAPA.
}
\usage{
mapafor(y, mapafit, fh = -1, ifh = 1, 
        comb = c("w.mean","w.median","mean","median","wght"), 
        outplot = c(0,1), hybrid = c(TRUE, FALSE), 
        conf.lvl = NULL, xreg=NULL)
}
\arguments{
  \item{y}{
    In sample observations of a time series (vector).
}
  \item{mapafit}{
    Fitted MAPA model (from mapaest).
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{ifh}{
    In-sample forecast horizon. Default = 0.
}
  \item{comb}{
    Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 1. 
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{conf.lvl}{
    Vector of confidence level for prediction intervals. Values must be (0,1). If conf.lvl == NULL then no intervals are calculated. For example to get the intervals for 80\% and 95\% use conf.lvl=c(0.8,0.95).
}
  \item{xreg}{
    Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample plus fh. Additional observations are unused. 
}
}
\value{
\item{infor}{In-sample forecasts.}
\item{outfor}{Out-of-sample forecasts.}
\item{PI}{Prediction intervals for given confidence levels.}
\item{MSE}{In-sample MSE error.}
\item{MAE}{In-sample MAE error.}
}
\note{The calculation of the prediction intervals is based on the empirical multiple step ahead MSE. To speed up calculations set conf.lvl=NULL. If very long forecast horizons are requested then once no more t+h MSE can be calculated the following approximation is used: sqrt(MSE_{t+1})*sqrt(h) for the error.}
\references{
    \itemize{
        \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time 
series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
        \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
        \item{You can find more information about MAPA at Nikos' \href{http://kourentzes.com/forecasting/tag/mapa/}{blog}.}
  }
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}; Fotios Petropoulos.
}
\seealso{
\code{\link{mapa}}, \code{\link{mapaest}}, \code{\link{mapacalc}}.
}
\examples{
mapafit <- mapaest(admissions,outplot=0)
out <- mapafor(admissions,mapafit)
}
\keyword{ ~mapa }
\keyword{ ~mapacalc }
\keyword{ ~mapaest }
