\name{tsaggr}
\alias{tsaggr}

\title{Non-overlapping temporal aggregation}

\description{Non-overlapping temporal aggregation.}

\usage{
tsaggr(y,fout,fmean=c(FALSE,TRUE),outplot=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    Time series vector (can be ts object).
}
  \item{fout}{
    Vector containing desirable aggregation levels. Must be positive and integer. If larger than length(y) then it is ignored.
}
  \item{fmean}{
    If TRUE the aggregated is done using mean, otherwise sum is used.
}
  \item{outplot}{
    If TRUE a plot of the original series and the aggregated ones is produced.
}
}
\value{
\item{out}{
  List of temporally aggregated series. If y was a ts object, then 'out' has ts objects with appropriate frequencies. Any non-integer frequency is set equal to 1. Series are named ALx, where x is the aggregation level. 
}
\item{all}{
  An array containing all aggregated series in the original frequency. Series are named ALx, where x is the aggregation level.
}
\item{idx}{
  List of indices used to produce 'out' from 'all': y.out[[i]] <- y.all[y.idx[[i]],i]. Series are named ALx, where x is the aggregation level.
}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
}
\examples{
out <- tsaggr(admissions,fout=2:12,fmean=TRUE,outplot=TRUE)
}
\keyword{ ~tsaggr }