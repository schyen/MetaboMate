#' Equidistant binning of spectra
#' @aliases binning
#' @export
#' @param X NMR matrix, rows represent spectra.
#' @param ppm Matching ppm vector.
#' @param width Width of bins
#' @param npoints Desired number of data points in binned ppm vector
#' @return Returned is a matrix with spectra in rows and binned ppm variables in columns. Column names represent binned ppm variables.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
binning=function(X, ppm, width=0.001, npoints=NULL){

  if(ppm[1]<ppm[length(ppm)]){stop('Let\'s stick to the convention: Revert the order of ppm and NMR matrix, so that the chemical shift decreases with increasing indices!')}
  if(ncol(X)!=length(ppm)){stop('Ppm vector does not match to NMR matrix dimension.')}

  if(is.null(width) & is.null(npoints)){
    stop('Pleas specify spectral bin width or the desired number of data points per spectrum.\n')
  }

  if(!is.null(width) & !is.null(npoints)){
    cat('Setting spectral width to 0.005 ppm.\n')
    width=0.005
  }

  if(!is.null(width)){
    if(width<=abs(diff(ppm[1:2]))){stop('Bin width equals or is smaller than the difference of neighbouring ppm points.')}
    ppm_bin=seq(min(ppm), max(ppm), by=width)
  }
  if(!is.null(npoints)){
    if(npoints>=length(ppm)){stop('Input variable npoints cannot be larger or equal than length of ppm vector.')}
    ppm_bin=seq(min(ppm), max(ppm), length.out=npoints)
  }

  seq(ppm[1], )
  quot=ceiling(width/(ppm[1]-ppm[2]))
  rr=ceiling(length(ppm)/quot)
  x=rep(1:200, each=rr)
  x=x[1:length(ppm)]



  Xb <- t(apply(X, 1, function(x, ppmt = ppm_bin, ppm_fres=ppm) {
    sInter = approxfun(ppm_fres ,x)
    sInter(ppmt)
  }))
  colnames(Xb) <- ppm_bin
  rownames(Xb) <- rownames(X)
  return(Xb)
}
