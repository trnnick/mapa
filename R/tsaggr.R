#' Non-overlapping temporal aggregation
#'
#' @description Non-overlapping temporal aggregation.
#'
#' @param y Time series vector (can be ts object).
#' @param fout Vector containing desirable aggregation levels. Must be positive and integer. If larger than length(y) then it is ignored.
#' @param fmean If TRUE the aggregated is done using mean, otherwise sum is used.
#' @param outplot If TRUE a plot of the original series and the aggregated ones is produced.
#' @param col.AL Vector of colours for the aggregated series in the plot. If NULL the spectral colouring scheme is used, with red the lowest and blue the highest aggregation level, as in the component plots of \code{\link{plot.mapa.calc}}.
#' @param ... Additional arguments passed to \code{plot}, for example \code{main}, \code{xlab}, and \code{ylab}.
#'
#' @return
#' \item{out}{
#'   List of temporally aggregated series. If y was a ts object, then 'out' has ts objects with appropriate frequencies. Any non-integer frequency is set equal to 1. Series are named ALx, where x is the aggregation level.
#' }
#' \item{all}{
#'   An array containing all aggregated series in the original frequency. Series are named ALx, where x is the aggregation level.
#' }
#' \item{idx}{
#'   List of indices used to produce 'out' from 'all': y.out[[i]] <- y.all[y.idx[[i]],i]. Series are named ALx, where x is the aggregation level.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@@kourentzes.com}.
#'
#' @keywords ~tsaggr
#' @examples
#' # Aggregate to quarterly and annual series, using sums
#' out <- tsaggr(admissions, fout=c(3,12))
#' out$out$AL12
#'
#' # All aggregation levels from 2 to 12, using means, with a plot
#' out <- tsaggr(admissions, fout=2:12, fmean=TRUE, outplot=TRUE, main="Admissions")
#'
#' # User defined colours for the aggregated series
#' out <- tsaggr(admissions, fout=c(3,12), fmean=TRUE, outplot=TRUE, col.AL=c("blue","red"))
#' @export
tsaggr <- function(y,fout,fmean=c(FALSE,TRUE),outplot=c(FALSE,TRUE),col.AL=NULL,...){
  # Defaults
  fmean <- fmean[1]
  outplot <- outplot[1]
  
  n <- length(y)
  
  # Get original frequency if y is ts
  if (isa(y,"ts")){
    f <- frequency(y)
  } else {
    f <- NULL
  }
  
  # Check fout
  fout <- fout[fout>0]
  fout <- fout[fout<=n]
  if (length(fout)==0){
    stop('fout must be positive integer(s), smaller than length(y).')
  }
  k <- length(fout)
  
  # Aggregation scale
  if (fmean == TRUE){
    scl <- fout
  } else {
    scl <- rep(1,k)
  }
  
  # Initialise array with aggregate series
  y.all <- array(NA,c(n,k))
  y.idx <- vector("list",k)
  y.out <- vector("list",k)
  slost <- vector("numeric",k)
  
  for (i in 1:k){
    # Find number of observations to drop from beginning
    slost[i] <- n %% fout[i]
    y.idx[[i]] <- seq((slost[i]+1),(n-fout[i]+1),fout[i])
    # Aggregate
    yk <- colSums(matrix(y[(slost[i]+1):n],nrow=fout[i]))/scl[i]
    y.all[(slost[i]+1):n,i] <- rep(yk,each=fout[i])
    # Produce aggregate series
    y.out[[i]] <- y.all[y.idx[[i]],i]
    if (!is.null(f)){
      if (f %% fout[i] == 0){
        ftemp <- f/fout[i]
      } else {
        ftemp <- 1
      }
      y.out[[i]] <- ts(y.out[[i]],frequency=ftemp)
    }
  }
  
  # Name outputs
  names(y.out) <- paste0("AL",fout)
  names(y.idx) <- paste0("AL",fout)
  colnames(y.all) <- paste0("AL",fout)
  
  # Produce plot
  if (outplot==TRUE){
    # Same colouring scheme as the MAPA component plots
    # cmp <- rainbow(k,start=0/6,end=4/6)
    if (is.null(col.AL)){
      cmp <- colorRampPalette(brewer.pal(11,"Spectral")[c(1:5,8:11)])(k)
    } else {
      cmp <- rep(col.AL,length.out=k)
    }
    # User arguments override defaults
    args <- list(x=as.numeric(y),col="black",type="l",ylab="",xlab="Period",lwd=2,
                 xaxs="i",yaxs="i",ylim=mapa.ylim(y,y.all))
    args <- modifyList(args,list(...))
    do.call(plot,args)
    for (i in 1:k){
      lines(y.all[,i],col=cmp[i],type="l")
      points((y.idx[[i]]+fout[i]-1),y.all[y.idx[[i]],i],pch=20,col=cmp[i])
    }
    legend("topleft",c("Original",paste0("AL",fout[c(1,k)])),col=c("black",cmp[1],cmp[k]),lty=1,lwd=c(2,1,1),cex=0.7)
  }

  return(list(out=y.out,all=y.all,idx=y.idx))
  
}