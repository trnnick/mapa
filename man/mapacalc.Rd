\name{mapacalc}
\alias{mapacalc}
\title{MAPA Internal Calculation of Forecasts}
\description{
The following function produces a simple MAPA forecast from a given origin. 
This is meant to be an internal function. Use mapafor instead.
}
\usage{
mapacalc(y, mapafit, fh = 0, comb = c("w.mean","w.median","mean","median","wght"), 
         outplot = c(0,1,2), hybrid = c(TRUE,FALSE), xreg=NULL)
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
  \item{comb}{
    Combination operator. This can be: "mean"; "median"; "wght" - where each aggregation level is weighted inversly to aggregation; "w.mean" - level and trend components are averaged, but seasonal and xreg follow the wght combination; "w.median" - as w.mean, but with median. It is suggested that for data with high sampling frequency to use one of the "w.mean" and "w.median".
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, forecasts and components. For the components the spectral colouring scheme is used. Dark red is aggregation level 1. Default is 1. 
}
\item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{xreg}{
    Vector or matrix of exogenous variables to be included in the MAPA. If matrix then rows are observations and columns are variables. Must be at least as long as in-sample plus fh. Additional observations are unused. 
}
}
\value{
\item{forecasts}{Vector with forecasts.}
\item{components}{Array with MAPA components.}
}
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
\code{\link{mapafor}}, \code{\link{mapa}}.
}
\examples{
mapafit <- mapaest(admissions,outplot=0)
mapacalc(admissions,mapafit,outplot=2)
}
\keyword{ ~mapafor }
