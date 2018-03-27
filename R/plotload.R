#' Plotting PCA or OPLS loadings
#' @export
#' @param model PCA or OPLS model generated via \emph{MetaboMate} package functions.
#' @param X Input matrix with rows representing spectra
#' @param ppm ppm variable
#' @param shift ppm region to visualise.
#' @param pc index of principal component to visualise, set to 1 if input model is OPLS
#' @param type Type of loadigs visualisation, either \code{'Statistical reconstruction'} or \code{'Backscaled'} (see Details).
#' @param title Plot title.
#' @details OPLS: If \code{type='Statistical recostruction'} the function calculates the covariance (y axis) and Pearson's correlation (colouring) of the predictive OPLS scores with each X variable (x axis is ppm varaible). If \code{type='Backscaled'} the OPLS loadings are backscaled with X feature standard deviations. Results are plotted over ppm, coloured according to OPLS model weights. Often, the latter method visualises model importance more robust due to the presence of false positive correlations. PCA: Function always calculates the statistical recostruction.
#' @author Torben Kimhofer
#' @importFrom stats cor cov
#' @importFrom ggplot2 ggplot geom_line scale_x_reverse scale_color_gradientn ggtitle xlab ylab theme_bw ggtitle aes_string
#' @importFrom colorRamps matlab.like2
#' @seealso \code{\link{OPLS_MetaboMate-class}}
#' @seealso \code{\link{opls}}
#' @seealso \code{\link{PCA_MetaboMate-class}}
#' @seealso \code{\link{pca}}


plotload=function(model, X, ppm, shift=c(0,10), pc=1, type=c('Statistical reconstruction', 'Backscaled'), title=''){

  if(length(model@Xscale)!=ncol(X)){
    stop('Model loadings do not fit to X matrix.')
  }

  if(length(model@Xscale)!=length(ppm)){
    stop('Model loadings do not fit to ppm vector.')
  }

  if(class(model)=='PLS_MataboMate'){
    type=c('Statistical reconstruction')
  }

  idx=get.idx(shift, ppm)

  if(type=='Statistical reconstruction'){
    switch(class(model)[1],
           "PCA_MetabMate"={t=model@t[,pc]},
           "OPLS_MetabMate"={t=model@t_pred[,pc]}
    )

    cc=cor(t, X)[1,]
    raCol=c(0, max(abs(cc)))
    cv=cov(t, X)[1,]

    df=data.frame(cor=abs(cc[idx]), cov=cv[idx], ppm=ppm[idx])

    g=ggplot(df, aes_string('ppm', 'cov', colour='cor'))+
      geom_line()+
      scale_x_reverse(breaks=seq(shift[1], shift[2], by=abs(diff(shift))/5))+
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
      scale_x_reverse(breaks=seq(shift[1], shift[2], by=abs(diff(shift))/5))+
      scale_color_gradientn(colors=matlab.like2(10),limits=raCol, name=expression(abs~w[pred*','~sc]))+
      ggtitle(title)+
      xlab(expression(delta~{}^1*H~(ppm)))+
      ylab(expression(p[pred]~sigma[x]))+
      theme_bw()

  }


  return(g)
}
