#' Simple plotting of multiple NMR spectra overlayed
#' @export
#' @param ppm ppm vector.
#' @param X NMR matrix with spectra represented in rows.
#' @param shift Chemical shift region to be plotted (in ppm).
#' @param add Logical indicating if spectra should be added to a current plot generated with \code{spec()} or \code{matspec()}.
#' @param ... Additional parameters to be passed on to the plot function.
#' @seealso \code{\link{spec}} \code{\link{plot}}
#' @aliases matspec
#' @details Low-level plotting function for NMR spectra.
#' @importFrom graphics matplot matpoints
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
matspec <- function(ppm, X, shift = c(0, 9.5), add = F, ...) {
    idx <- get.idx(shift, ppm)
    if (add == T) {
        matpoints(ppm[idx], t(X[, idx]), type = "l", ...)
    } else {
        matplot(ppm[idx], t(X[, idx]), type = "l", xlim = rev(range(ppm[idx])), xlab = "ppm", ylab = "Intensity (AU)", ...)
    }
}
