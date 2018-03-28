#' Probablistic quotient normalisation
#' @export
#' @aliases pqn
#' @param X Input matrix where rows represent samples.
#' @param reference.idx Indices of spectra in X used to calculate reference spectrum (see Details).
#' @param TArea Logical indicating if total area normalisation should be applied first (see Details).
#' @param add.DilF Character string for new variable for dilution factor (will be exported as global envrionment variable). Can be left NULL if dilution factor should not be exported.
#' @details It is sometimes favourable not to use all spectra to calculate a dilution reference (e.g. QC samples should generally be excluded). Therefore, a vector of indices can be specified with the parameter \code{reference.idx} and the respective spectra are used to calculate the median spectrum as a dilution reference. The parameter \code{reference.idx} can also be a single index, then the respective spectrum is used as a reference. If it is set \code{N/A}, all spectra in \code{X} are used to calculate the dilution reference spectrum. to Total area normaliasation is integral part of the probablistic quotient normalisation algorithm (see References), however, this sometimes distorts the spectra, i.e. is not suitable for normalisation. The parameter \code{TArea} can be set to TRUE or FALSE, depending if total area normalisation should be applied or not.
#' @references Dieterly, F., \emph{et al.} (2006), Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. Application in 1H NMR Metabonomics, \emph{Analytical Chemistry}, 78.3, 4281-90.
#' @author Torben Kimhofer
pqn=function(X, reference.idx=NA, TArea=F, add.DilF=NULL){

  # apply total area normalisation
  if(TArea==T){
    X=t(apply(X, 1, function(x) x/sum(x)*1000))
  }

  le.ref=length(reference.idx)

  if(le.ref==1 & !is.na(reference.idx[1])) {ref=X[reference.idx,]}
  if(le.ref==1 & is.na(reference.idx[1])) {ref=apply(X, 2, median)}
  if(le.ref>1) {ref=apply(X[reference.idx,], 2, median)}

  ref[ref==0]=0.0001

  dil.F=1/apply(X, 1, function(x) median(x/ref))
  X.pqn=t(apply(rbind(1:nrow(X)), 2, function(i){
    X[i,]*dil.F[i]
  }))

  if(!is.null(add.DilF)){
    assign(add.DilF, dil.F, envir=.GlobalEnv)
  }

  return(X.pqn)
}
