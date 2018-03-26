#' Spectral calibration to a chemical shift reference
#' @export
#' @aliases calibration
#' @description Spectral calibration to a chemical shift reference.
#' @param NMR Input NMR matrix with rows representing spectra and columns representing chemical shifts.
#' @param ppm ppm vector with the same length as \code{ncol(X)}.
#' @param type Either 'Urine' or 'Plasma' for urine or blood-derived spectra, respectively (see Details).
#' @details Spectral calibration to a chemical shift reference. If \code{type='Urine'} calibration will be performed using signal resonating around 0 ppm (Trimethylsilylpropanoic acid resonance). If \code{type='Plasma'}, the doublet around 5.23, originating from the alpha anomer of glucose, will be used for calibration. \strong{Blood serum-derived spectra} can also be calibrated with input argument \code{type='Plasma'}.
#' @return Returned is the calibrated NMR data matrix.
#' @references Dona, A.C., \emph{et al.} (2014) Precision high-throughput proton NMR spectroscopy of human urine, serum, and plasma for large-scale metabolic phenotyping. \emph{Analytical Chemistry}. 86.19. 9887-94.
#' @author Torben Kimhofer
#' @importFrom  speaq detectSpecPeaks

### function shiftes the max TSP intensity to zero
calibration=function(NMR, ppm, type='Urine'){
  # function get indices in ppm which match the specified range
  get.idx=function(range=c(1,5), ppm){
    range=sort(range, decreasing=T);
    which(ppm<=range[1] & ppm>=range[2])}

  if(class(type[1])=='numeric'){
    idx=get.idx(c(type[1:2]), ppm)
    zero.ppm=which.min(abs(ppm[idx]-type[3]))

    maxInt=array()
    for(i in 1:nrow(NMR)){
      maxInt[i]=which.max(NMR[i,idx])
    }

  }else{

  if(type=='Urine'){
    idx=get.idx(c(-0.2, 0.2), ppm)
    zero.ppm=which.min(abs(ppm[idx]))

    maxInt=array()
    for(i in 1:nrow(NMR)){
      maxInt[i]=which.max(NMR[i,idx])
    }

  }

  if(type=='Plasma'){
    idx=get.idx(c(5.1, 5.3), ppm)
    zero.ppm=which.min(abs(5.233-(ppm[idx])))

    # perform peak picking (squaq, inlcluding noise estimation)
    # find two signals that have a certain distance and are close to 5.233
    source('/Users/tk2812/Box Sync/R_scr/fct_pqn.R')
    tt=noise.est(NMR, ppm, where=c(14.4, 14.5))

    ppm1=ppm[idx]
    X1=NMR[,idx]
    res=list()
    for(i in 1:nrow(X1)){
      res[[i]]=detectSpecPeaks(t(matrix(X1[i,])), baselineThresh=tt[i]*4)
    }


    ics=list()
    for(i in 1:nrow(X1)){
      idx=which(diff(res[[i]][[1]])>39 & diff(res[[i]][[1]])<44)
      if(length(idx)>1){
        distGluc=array()
        for(j in 1:length(idx)){
          distGluc[j]=mean(ppm1[res[[i]][[1]][idx[j]]:res[[i]][[1]][idx[j]+1]])
        }
        id=which.min(abs(distGluc-5.233))
        ics[[i]]=res[[i]][[1]][idx[id]:(idx[id]+1)]
      }else{
        ics[[i]]=res[[i]][[1]][idx:(idx+1)]
      }
    }
    if(length(which(sapply(ics, length)>2))>0){stop('Something went wrong!')}
    maxInt=sapply(ics, '[[', 1)
  }
  }

  Int.corr=zero.ppm-maxInt
  # if Int.corr<0: shift is > zero

  for(i in 1:nrow(NMR)){

   x=NMR[i,]

   if(Int.corr[i]<0){
     x=c(x[-c(1:abs(Int.corr[i]))], rep(0, abs(Int.corr[i])))
   }
   if(Int.corr[i]>0){
     x=c(rep(0, abs(Int.corr[i])),x)
   }

   NMR[i,]=x[1:length(ppm)]

  }

  return((NMR))
}
