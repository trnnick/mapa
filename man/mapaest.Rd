\name{mapaest}
\alias{mapaest}
\title{MAPA Estimation}
\description{
The following function estimates MAPA and saves all fitted model.
}
\usage{
mapaest(y, ppy, minimumAL = 1, maximumAL = ppy, paral = c(0,1,2), 
	      display = c(0,1), outplot = c(0,1), model = "ZZZ", type=c("ets","es"), ...)
}
\arguments{
  \item{y}{
    In sample observations of a time series (vector).
  }
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
  }
  \item{minimumAL}{
    Lowest aggregation level to use. Default = 1, maximumAL>1.
  }
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy.
  }
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes 
	and initialise cluster. Default is 0.
  }
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
  }
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 0.
  }
  \item{model}{
    Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M" or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md", "Z", "X" or "Y"); and the third letter refers to the season type ("N","A","M", "Z", "X" or "Y"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative, "Z"=automatically selected, "X"=automatically select between none and additive and "Y"=automatically select between none and multiplicative. A "d" for trend implies damped. "X" and "Y" supported only by type=="es". If used with type=="ets" a warning will be given and they will default to "Z". By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. For aggregation levels that seasonality becomes 1 then a non-seasonal model is estimated.
  }
  \item{type}{
    What type of exponential smoothing implementation to use. "es" = use from the smooth package; "ets" = use from the forecast package. Default is "es"
  }
  \item{...}{
    Pass additional arguments to es or ets.
  }
}
\value{
\item{mapafit}{Estimated MAPA model structure.}
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
\code{\link{mapa}}, \code{\link{mapafor}}.
}
\examples{
mapafit <- mapaest(admissions)
}
\keyword{ ~mapa }
\keyword{ ~mapafor }