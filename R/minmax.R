#' Min-Max scaling
#' @export
#' @details Scaling minimum and maximum values to zero and one
#' @param x Input vector.
#' @return Scaled x vector.
#' @usage minmax(x)
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
minmax<-function(x){
  (x-min(x))/(max(x)-min(x))}
