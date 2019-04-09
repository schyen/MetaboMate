#' Plotting PCA or OPLS loadings
#' @export
#' @param model PCA or OPLS model generated via \emph{MetaboMate} package functions.
#' @param X Input matrix with rows representing spectra
#' @param ppm ppm variable
#' @param shift ppm region to visualise.
#' @param pc index of principal component to visualise, set to 1 if input model is OPLS
#' @param type Type of loadings visualisation, either \code{'Statistical reconstruction'} or \code{'Backscaled'} (see Details).
#' @param title Plot title.
#' @details OPLS: If \code{type='Statistical reconstruction'} the function calculates the covariance (y axis) and Pearson's correlation (colouring) of the predictive OPLS scores with each X variable (x axis is ppm variable). If \code{type='Backscaled'} the OPLS loadings are backscaled with X feature standard deviations. Results are plotted over ppm, coloured according to OPLS model weights. Often, the latter method visualises model importance more robust due to the presence of false positive correlations. PCA: Function always calculates the statistical recostruction.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @references Cloarec, O., \emph{et al.} (2005). Evaluation of the Orthogonal Projection on Latent Structure Model Limitations Caused by Chemical Shift Variability and Improved Visualization of Biomarker Changes in 1H NMR Spectroscopic Metabonomic Studies. \emph{Analytical Chemistry} 77.2, 517-26.
#' @importFrom stats cor cov
#' @importFrom ggplot2 ggplot geom_line scale_x_reverse scale_color_gradientn ggtitle xlab ylab theme_bw ggtitle aes_string
#' @importFrom colorRamps matlab.like2
#' @importFrom scales pretty_breaks
#' @seealso \code{\link{pca}} \code{\link{opls}} \code{\link{PCA_MetaboMate-class}} \code{\link{OPLS_MetaboMate-class}}

plotload=function(model, X, ppm, shift=c(0,10), pc=1, type=c('Statistical reconstruction', 'Backscaled'), title=''){

  if(grepl('stat|recon', type, ignore.case = T)){type='Statistical reconstruction'}else{
    type='Backscaled'
  }

  # # why is this in here?
  if(class(model)[1]=='PCA_MetaboMate'){
    #type=c('Statistical reconstruction')


    if(nrow(model@p)!=ncol(X)){
      stop('Model loadings do not fit to X matrix.')
    }

    if(nrow(model@p)!=length(ppm)){
      stop('Model loadings do not fit to ppm vector.')
    }

  }else{

    if(ncol(model@p_pred)!=ncol(X)){
      stop('Model loadings do not fit to X matrix.')
    }

    if(ncol(model@p_pred)!=length(ppm)){
      stop('Model loadings do not fit to ppm vector.')
    }

  }

  idx=get.idx(shift, ppm)

  if(type=='Statistical reconstruction'){
    switch(class(model)[1],
           "PCA_MetaboMate"={t=model@t[,pc]},
           "OPLS_MetaboMate"={t=model@t_pred[,pc]}
    )

    cc=cor(t, X)[1,]
    raCol=c(0, max(abs(cc)))
    cv=cov(t, X)[1,]

    df=data.frame(cor=abs(cc[idx]), cov=cv[idx], ppm=ppm[idx])

    g=ggplot(df, aes_string('ppm', 'cov', colour='cor'))+
      geom_line()+
      scale_x_reverse(breaks=pretty_breaks(n=15))+
      scale_color_gradientn(colors=matlab.like2(10), limits=raCol, name='cor(t,x)')+
      ggtitle(title)+
      xlab(expression(delta~{}^1*H~(ppm)))+
      ylab('cov(t,x)')+
      theme_bw()
  }

  if(type=='Backscaled'){

    # backscaling p
    p=model@p_pred[pc,]*model@Xscale
    w=minmax(abs(model@w_pred[pc,]))

    raCol=range(w)
    df=data.frame(t_bs=p, w=abs(w), ppm)
    df=df[idx,]

    g=ggplot(df, aes_string('ppm', 't_bs', colour='w'))+
      geom_line()+
      scale_x_reverse(breaks=pretty_breaks(n=15))+
      scale_color_gradientn(colors=matlab.like2(10),limits=raCol, name=expression(abs~w[pred*','~sc]))+
      ggtitle(title)+
      xlab(expression(delta~{}^1*H~(ppm)))+
      ylab(expression(p[pred]~sigma[x]))+
      theme_bw()

  }


  return(g)
}
