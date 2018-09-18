#' get indices of respective ppm values
#' @param range ppm range
#' @param ppm ppm vector
#' @export
#' @aliases get.idx
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
get.idx=function(range=c(1,5), ppm){
  range<-sort(range, decreasing=T);
  which(ppm<=range[1] & ppm>=range[2])
  }
