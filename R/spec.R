#' Plotting a single NMR spectrum
#' @export
#' @param ppm ppm vector.
#' @param x NMR spectrum.
#' @param shift chemical shift region to be plotted.
#' @param add Logical indicating if spectrum should be added to a current plot generated with \code{spec()} or \code{matspec()}.
#' @param ... Additional parameters to be passed on to the plot function.
#' @seealso  \code{\link{matspec}} \code{\link{plot}}
#' @aliases spec
#' @details Low-level plotting function for a single NMR spectrum.
#' @importFrom graphics points

spec=function(ppm, x, shift=c(0,9.5), add=F, ...){
  idx=get.idx(shift, ppm)
  if(add==T){
    points(ppm[idx], x[idx], type='l',  ...)
  }else{
    plot(ppm[idx], x[idx], type='l', xlim=rev(range(ppm[idx])), xlab='ppm', ylab='Intensity (AU)', ...)}

}
