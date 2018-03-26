#' Centering and scaleing vectors or matrices
#' @param X Data to be scaled, can be a vector or a matrix with observations in rows variables in columns.
#' @param idc Indices of samples in X used to derive center/scale parameters (relevant for statistical validation procedures, e.g., scaling X with training set samples). If set to \code{all}, no selection is performed.
#' @param center Logical indicating if data should be mean centered.
#' @param scale Logical indicating if data should be scaled. Currently only unit variance (UV) scaling is implemented.
#' @return Centred/scaled data of the same dimensions than input.


center_scale <- function(X, idc='all', center=T, scale = 'UV') {

  # if idc empty then use all samples
  if(idc[1]=='all'){
    X1<-X }else{
    if(ncol(X)==1){
      X1<-cbind(X[idc,]) }else{
        X1<-X[idc,] }
  }

  # centering and scaling, use parameters from training set (idc)
  if (is.logical(center[1]) & center[1]==T ) {
    meanX <- apply(X1, 2, mean)
    X <- scale(X, center = meanX, scale = F)
  }

  if (is.numeric(center[1])) {
    X <- scale(X, center = center, scale = F)
  }

  if (scale[1] == 'UV') {
    sdX <- apply(X1, 2, sd)
    X <- scale(X, center = F, scale = sdX)
  }

  if (is.numeric(scale[1])) {
    X <- scale(X, center = F, scale = scale)
  }

  return(X)
}
