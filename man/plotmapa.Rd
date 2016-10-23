\name{plotmapa}
\alias{plotmapa}
\title{Produce estimated MAPA fit plot (OBSOLETE!)}
\description{
The following function plots the fitted components at each aggregation
level of an estimated MAPA. The plot is the same as the one that can be 
produced optionally by mapaest.
}
\usage{
plotmapa(mapafit)
}
\arguments{
  \item{mapafit}{
    Fitted MAPA model (from mapaest).
}
}
\references{
Kourentzes N., Petropoulos F., Trapero J.R., 2014. Improving forecasting by estimating time 
series structural components across multiple frequencies, International Journal of Forecasting,
 30(2), 291-302.
}
\author{
Nikolaos Kourentzes
}
\seealso{
\code{\link{mapaest}}, \code{\link{mapa}}.
}
\examples{
mapafit <- mapaest(admissions,outplot=0)
plotmapa(mapafit)
}
\keyword{ ~plotmapa }
