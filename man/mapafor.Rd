\name{mapafor}
\alias{mapafor}
\title{MAPA Forecast}
\description{
The following function produces in- and out-of-sample MAPA forecasts, for multiple steps ahead.
 This is the recommended function to use in forecasting with MAPA.
}
\usage{
mapafor(y, mapafit, fh = -1, ifh = 1, comb = c("mean","median","wght"), 
        outplot = c(1,0), hybrid = c(TRUE, FALSE), 
        conf.lvl = NULL)
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
    Combination operator. One of "mean" or "median". Default is "mean". "wght" is experimental, use at your own risk!
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
Kourentzes N., Petropoulos F., Trapero J.R., 2014. Improving forecasting by estimating time 
series structural components across multiple frequencies, International Journal of Forecasting,
 30(2), 291-302.
}
\author{
Nikolaos Kourentzes and Fotios Petropoulos
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