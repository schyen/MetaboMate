#' Calculating full width at half max
#' @export
#' @param X NMR matrix with rows representing spectra.
#' @param ppm ppm vector with its length equals to \code{nrow(X)}.
#' @param shift Signal shift to calculate line width
#' @description Calculating full width at half maximum (FWHM, aka line width). This function simply returns the ppm difference where peak line crosses half of the peak height. It requires one signal across all spectra within ppm ranges specified in \code{shift}.
#' @return Array of line widths in ppm. To convert from ppm to Hertz (Hz), multiply values with the spectrometer frequency (column \code{a_SF01} in \code{meta} data frame).
#' @seealso \code{\link{readBruker}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
lw <- function(X, ppm, shift = c(-0.01, 0.01)) {
    idx <- get.idx(shift, ppm)
    fwhm <- apply(X[, idx], 1, function(x, pp = ppm[idx]) {
        x <- x + abs(min(x))  # no negative values
        baseline <- min(x)
        height <- max(x) - baseline
        hw <- baseline + (0.5 * height)
        f <- approxfun(pp, x)
        x_new <- seq(-0.01, 0.01, by = 1e-05)
        y_new <- f(x_new)
        diff(x_new[range(which(y_new > hw))])
    })
    return(fwhm)
}
