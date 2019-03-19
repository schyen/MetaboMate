#' Min-max scaling
#' @export
#' @param x Numeric vector to be scaled.
#' @return Scaled x vector.
#'
#' @details Data is scaled to range between zero and one:
#' \deqn{X_{scaled}=\frac{x-x_{min}}{x_{max}-x_{min}}}
#' @usage minmax(x)
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
minmax <- function(x) {
    (x - min(x))/(max(x) - min(x))
}
