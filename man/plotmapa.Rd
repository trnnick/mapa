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
\note{The function 'plotmapa' is deprecated. Use plot() instead.}
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
\code{\link{mapaest}}, \code{\link{mapa}}.
}
\examples{
mapafit <- mapaest(admissions,outplot=0)
plotmapa(mapafit)
}
\keyword{ ~plotmapa }
