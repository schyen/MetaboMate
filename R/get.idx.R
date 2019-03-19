#' Find indices in ppm vector for respective chemical shift range
#' @param range Chemical shift range (in ppm)
#' @param ppm Ppm vector
#' @export
#' @aliases get.idx
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
get.idx <- function(range = c(1, 5), ppm) {
    range <- sort(range, decreasing = T)
    which(ppm <= range[1] & ppm >= range[2])
}
