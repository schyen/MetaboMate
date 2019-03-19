#' Estimation of noise level
#' @export
#' @param NMR Input matrix with rows representing spectra.
#' @param ppm ppm vector with its length equals to \code{ncol(X)}.
#' @param where Signal free region across all NMR spectra (see Details).
#' @details Estimation of noise level in NMR spectra. This is useful for quality control checks (e.g., before and after spectral normalisation). Noise estimation requires a signal-free ppm region across all spectra, usually this is at the extreme ends or the spectrum. This function requires a minimum number of 50 data points to robustly estimate noise levels.
#' @return Returned is a vector of noise levels for each spectrum.
#' @seealso \code{\link{readBruker}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom ptw asysm
noise.est <- function(NMR, ppm, where = c(14.6, 14.7)) {
    # set ppm range where noise should be estimated, i.e., no signals
    idx <- get.idx(where, ppm)
    if (length(idx) < 50) {
        stop("No or too few data points for noise estimation!")
    }
    noise <- apply(NMR[, idx], 1, function(x) {
        x_driftcorrected <- x - asysm(x, lambda = 1e+10)
        noi <- (max(x_driftcorrected) - min(x_driftcorrected))
        return(noi)
    })
    return(noise)
}
