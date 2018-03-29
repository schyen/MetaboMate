#' Overlay PCA or OPLS loadings with spectra
#' @export
#' @param model PCA or OPLS model generated via \emph{MetaboMate} package functions.
#' @param X Input matrix with rows representing spectra
#' @param ppm ppm variable
#' @param shift ppm region to visualise.
#' @param pc index of principal component to visualise, set to 1 if input model is OPLS
#' @param type Type of loadigs visualisation, either \code{'Statistical reconstruction'} or \code{'Backscaled'} (see Details).
#' @param an List with one to three elements specifying facetting, colour and linetype (see Details).
#' @param alp Alpha value for spectral lines.
#' @param title Plot title.
#' @param size plot line width.
#' @param ... Additional paramters passe on to ggplot's facet function.
#' @description  Plotting overlayed NMR specra. This function is based on ggplot2, a high-level plotting R package. For high ppm ranges computation time is relatively, so the range of input argument \code{shift} should be as small as possible. List argument \code{an} must have the first element define, even if it is only a single value. If colour and line width is specified, then at least one list elements of \code{an} must have the same length as \code{X}.
#' @details OPLS: If \code{type='Statistical recostruction'} the function calculates the covariance (y axis) and Pearson's correlation (colouring) of the predictive OPLS scores with each X variable (x axis is ppm varaible). If \code{type='Backscaled'} the OPLS loadings are backscaled with X feature standard deviations. Results are plotted over ppm, coloured according to OPLS model weights. Often, the latter method visualises model importance more robust due to the presence of false positive correlations. PCA: Function always calculates the statistical recostruction.
#' @seealso \code{\link{plotload}} \code{\link{specOverlay}} \code{\link[=OPLS_MetaboMate-class]{OPLS_MetaboMate}} \code{\link{opls}} \code{\link[=PCA_MetaboMate-class]{PCA_MetaboMate}} \code{\link{pca}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes_string scale_x_reverse ggtitle xlab facet_grid theme_bw theme element_text geom_line scale_colour_gradientn
#' @importFrom colorRamps matlab.like2
#' @importFrom scales pretty_breaks
#' @importFrom stats as.formula


specload=function(model, X, ppm, shift=c(0,10), an, alp=0.3, size=0.5, pc=1, type=c('Statistical reconstruction', 'Backscaled'), title=''){

  if(length(model@Xscale)!=ncol(X)){
    stop('Model loadings do not fit to X matrix.')
  }

  if(length(model@Xscale)!=length(ppm)){
    stop('Model loadings do not fit to ppm vector.')
  }

  if(class(model)=='PLS_MataboMate'){
    type=c('Statistical reconstruction')
  }

  if(is.numeric(an[[1]])){
    stop('Facet level is numeric!')
  }

  if(length(an)<3){
    if(length(an)==1){
      an[[2]]='black'
      an[[3]]='1'
    }else{
      an[[3]]='1'
    }
  }
  names(an)=c('facet', 'col', 'ltype')

  idx=get.idx(shift, ppm)
  if(length(idx)==0){stop('Shift area not in ppm variable.')}

  if(grepl('Stat|struct|ical', type, ignore.case = T)){
    type='Statistical reconstruction'
  }else{
    type='Backscaled'
  }

  ########## define overalay
  if(type=='Statistical reconstruction'){
    switch(class(model)[1],
           "PCA_MetaboMate"={t=model@t[,pc]},
           "OPLS_MetaboMate"={t=model@t_pred[,pc]}
    )

    cols=cor(t, X)[1,]
    raCol=c(0, max(abs(cols)))
    y=cov(t, X)[1,]

  }

  if(type=='Backscaled'){

    # backscaling p
    y=model@p_pred[pc,]*model@Xscale
    cols=minmax(abs(model@w_pred[pc,]))
    raCol=c(0,1)
  }

  #####################
  le.arg<-length(an)
  idx=get.idx(shift, ppm)

  specs=X[,idx]
  limY=range(specs)
  colnames(specs)=paste("ppm", ppm[idx], sep='_')

  df=data.frame(ID=1:nrow(specs), do.call(cbind.data.frame, an) , specs)

  df=melt(df, id.vars=c('ID', names(an)))
  df$variable=as.numeric(gsub('^\\.', '-', gsub('ppm_','', df$variable)))

  #colnames(df)[4:6]=c('Group', 'ppm', 'Intensity')
  df$load=NA


  # reduce to idx
  cols=abs(cols[idx])
  y=y[idx]
  cv1=(minmax(y)*(limY[2]/3))+limY[2]*0.67



  if(max(cv1)>limY[2]){
    cv1=cv1-abs(max(cv1-limY[2]))
  }


  fac.lev=unique(an[[1]])
  # define loadings
  df1=data.frame(alp,
                 ID=nrow(X)+1,
                 facet=fac.lev[1],
                 Group='load',
                 ppm=ppm[idx],
                 Intensity=cv1,
                 load=cols)

  for(i in 2:length(fac.lev)){
    df1=rbind(df1,
              data.frame(alp,
                         ID=nrow(X)+i,
                         facet=fac.lev[i],
                         Group='load',
                         ppm=ppm[idx],
                         Intensity=cv1,
                         load=cols))
  }

  g=ggplot()+
    geom_line(data=df1, aes_string('ppm', 'Intensity', color='load', group='ID'), size=0.8)+
    geom_line(data=df, aes_string('variable', 'value',group='ID'), alpha=alp, size=0.1)+
    scale_x_reverse(breaks=round(seq(shift[1], shift[2], by=abs(diff(shift))/20),3))+
    scale_y_continuous(limits=limY)+
    ggtitle(title)+
    xlab(expression(delta~{}^1*H~'(ppm)'))+
    ylab('Intensity (AU)')+
    facet_grid(facet~.)+
    theme_bw()+
    theme(axis.text = element_text(colour="black"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(type=='Statistical reconstruction'){
    g=g+scale_colour_gradientn(colors=matlab.like2(10), name='cor(t,x)', limits=raCol)
  }else{
    g=g+scale_colour_gradientn(colors=matlab.like2(10), name=expression(abs~w[pred*','~sc]), limits=raCol)
  }

  return(g)


}
