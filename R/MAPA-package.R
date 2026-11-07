# Package level documentation and imports for MAPA
# Nikolaos Kourentzes 2026
#
# The description repeats the Description field of DESCRIPTION, so keep the two 
# in sync. It is given here because roxygen converts the DOI in DESCRIPTION to 
# \doi{}, which is evaluated only when the package is built. Help previews from 
# the source files (for example devtools/RStudio) then show the raw \Sexpr code. 
# The \Sexpr below gives the same DOI link, but is evaluated when the help is 
# displayed, so it renders correctly everywhere.

#' @description Functions and wrappers for using the Multiple Aggregation Prediction Algorithm (MAPA) for time series forecasting. MAPA models and forecasts time series at multiple temporal aggregation levels, thus strengthening and attenuating the various time series components for better holistic estimation of its structure. For details see Kourentzes et al. (2014) \Sexpr[results=rd,stage=render]{tools:::Rd_expr_doi("10.1016/j.ijforecast.2013.09.006")}.
#'
#' @section Wrapper functions:
#' These estimate MAPA and produce forecasts in a single call.
#' \describe{
#'   \item{\code{\link{mapa}}}{Estimate MAPA and produce in- and out-of-sample forecasts, with optional prediction intervals. This is the recommended function for most uses.}
#'   \item{\code{\link{mapasimple}}}{Estimate MAPA and produce out-of-sample forecasts, together with the forecasted components at each aggregation level.}
#' }
#'
#' @section Estimation and forecasting:
#' These separate estimation from forecasting, so that a fitted MAPA can be reused, for example to forecast with new xreg values.
#' \describe{
#'   \item{\code{\link{mapaest}}}{Estimate the exponential smoothing model at each aggregation level, using \code{ets} from the \pkg{forecast} package or \code{es} from the \pkg{smooth} package, optionally with exogenous variables (MAPAx).}
#'   \item{\code{\link{mapafor}}}{Produce in- and out-of-sample forecasts from a fitted MAPA, with optional prediction intervals and in-sample errors.}
#'   \item{\code{\link{mapacalc}}}{Produce out-of-sample forecasts from a fitted MAPA, together with the forecasted components at each aggregation level.}
#' }
#'
#' @section Methods for MAPA outputs:
#' The outputs of the functions above have the classes \code{mapa.fit} (\code{mapaest}), \code{mapa.calc} (\code{mapacalc} and \code{mapasimple}), and \code{mapa.for} (\code{mapafor} and \code{mapa}), with the following methods.
#' \describe{
#'   \item{\code{\link{plot.mapa.fit}}}{Plot the components of the model fitted at each aggregation level.}
#'   \item{\code{\link{plot.mapa.calc}}}{Plot the forecasts and, optionally, the forecasted components at each aggregation level.}
#'   \item{\code{\link{plot.mapa.for}}}{Plot the in- and out-of-sample forecasts and the prediction intervals.}
#'   \item{\code{\link[=summary.mapa.fit]{print}} and \code{\link{summary.mapa.fit}}}{Print and summarise each of the outputs: the fitted model and parameters at each aggregation level; the combined components for each horizon; and the in-sample errors for each horizon.}
#' }
#'
#' @section Utilities:
#' \describe{
#'   \item{\code{\link{tsaggr}}}{Non-overlapping temporal aggregation of a time series to one or more aggregation levels.}
#' }
#'
#' @section Data:
#' \describe{
#'   \item{\code{\link{admissions}}}{Monthly UK NHS total non-elective G&A admissions, April 2008 to June 2013.}
#' }
#'
#' @section Compatibility with earlier versions:
#' These are provided only for compatibility and will be removed in a later version of the package.
#' \describe{
#'   \item{\code{\link{mapafit.old}}}{Transform a fitted MAPA to the structure of versions <= 2.0.6, for code that depends on the earlier indexing.}
#'   \item{\code{\link{plotmapa}}}{Plot a fitted MAPA. This is deprecated, use \code{plot} instead.}
#' }
#'
#' @references
#'     \itemize{
#'         \item{Kourentzes N., Petropoulos F., Trapero J.R. (2014) Improving forecasting by estimating time
#' series structural components across multiple frequencies. \emph{International Journal of Forecasting}, \bold{30}(\bold{2}), 291--302.}
#'         \item{Kourentzes N., Petropoulos F. (2015) Forecasting with multivariate temporal aggregation: The case of promotional modelling. \emph{International Journal of Production Economics}.}
#'         \item{You can find more information about MAPA at Nikos' \href{https://kourentzes.com/forecasting/tag/mapa/}{blog}.}
#'   }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}; Fotios Petropoulos.
#'
#' @keywords internal
#' @import forecast parallel RColorBrewer
#' @importFrom grDevices col2rgb rgb colorRampPalette
#' @importFrom graphics axis image layout legend lines par plot points polygon text box
#' @importFrom stats frequency median qnorm ts prcomp sd time tsp
#' @importFrom utils tail modifyList
#' @importFrom smooth es modelType
"_PACKAGE"
