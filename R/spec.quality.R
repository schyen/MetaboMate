#' Calculating spectral quality indices
#' @aliases spec.quality
#' @export
#' @param X NMR matrix, rows represent spectra.
#' @param ppm Matching ppm vector.
#' @param ppm.noise region for noise estimation, must be signal free across all spectra.
#' @param plot Logical indicating if a summary plot should be produced.
#' @return Returned is a data frame with TSP line width, residual water signal, normalised baseline estimation (the lower the less baseline fluctuations), estimated signal to noise ratio.
#' @references Eilers, P.H.C. (2004), Parametric Time Warping, \emph{Analytical Chemistry}, 76.2, 404–411.
#' @references Bloemberg, T.G., \emph{et al.} (2010), Improved parametric time warping for Proteomics, \emph{Chemometrics and Intelligent Laboratory Systems}, 104.1, 65-74.
#' @importFrom ptw asysm
#' @importFrom ggplot2 aes_string
#' @importFrom colorRamps matlab.like2
#' @importFrom graphics plot
#' @importFrom stats median
#'


spec.quality=function(X, ppm, ppm.noise=c(9.4,9.5), plot=T){

  # TSP line width estimataion
  tsp.lw=lw(X, ppm, shift=c(-0.01,0.01))

  # residual water
  idx=get.idx(range=c(4.7, 4.95), ppm)
  resW<-apply(X[,idx], 1, function(x, pp=ppm[idx]){
   sum(abs(x))
  })/length(idx)

  # baseline assessment
  idx=c(get.idx(range=c(0.1, 4.3), ppm), get.idx(range=c(5.2, 9.5), ppm))
  bl=apply(X[,idx], 1, function(x) {
    bl=asysm(x)
    #bl=bl+abs(min(bl))
    bl})

  bl.sum=apply(bl, 2, sum)/ncol(X)

  # noise estimation
  X.bl=t(apply(cbind(1:ncol(bl)), 1, function(i, idc=idx) {
    x.bl=(X[i,idc]-bl[,i])
    x.bl
    x.bl+abs(min(x.bl))
  }))
  spec.n=noise.est(X.bl, ppm[idx], where=c(9.4,9.5))

  # estimate signal to noise ratio
  sn.ratio=apply(cbind(1:ncol(bl)), 1, function(i){
    idx=X.bl[i,]>spec.n[i]
    median(X.bl[i,idx])/spec.n[i]
  })

  out=data.frame(TSP.lw.ppm=tsp.lw, Residual.water=resW, Baseline=bl.sum, SN.ratio=sn.ratio)

  if(plot==T){
    g=ggplot(data=out,
             aes_string(x='SN.ratio',
                        y='Residual.water',
                        colour='Baseline',
                        size='TSP.lw.ppm'))+
      geom_point()+
      scale_colour_gradientn(colours=matlab.like2(10))
    plot(g)
  }

  return(out)
}
