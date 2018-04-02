#' Overlay PCA or OPLS loadings with spectra - experimental
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
#' @importFrom ggplot2 aes_string scale_x_reverse ggtitle xlab facet_grid theme_bw theme element_text geom_line scale_colour_gradientn element_blank scale_colour_brewer
#' @importFrom colorRamps matlab.like2
#' @importFrom scales pretty_breaks
#' @importFrom stats as.formula
#' @importFrom grid grid.newpage viewport grid.layout pushViewport unit unit.length

specload1=function(model, X, ppm, shift=c(0,10), an, alp=0.3, size=0.5, pc=1, type=c('Statistical reconstruction', 'Backscaled'), title=''){


  if(class(model)=='PCA_MetaboMate'){
    type=c('Statistical reconstruction')

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

  if(is.numeric(an[[1]])){
    stop('Facet level is numeric!')
  }

  names(an)[1]='facet'

  if(is.null(names(an))){cat('No facet, colour and linetype names given. See an argument in ?specOverlay\n')
    names(an)=paste('an', 1:length(an), sep='')}

  if ('' %in% names(an)){
    idx=which(names(an)=='')
    names(an)[idx]=paste('an', idx, sep='')
  }
  names(an)=gsub(' ', '.', names(an))


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
                 size,
                 ID=nrow(X)+1,
                 facet=fac.lev[1],
                 Group='load',
                 ppm=ppm[idx],
                 Intensity=cv1,
                 load=cols)

  if(length(fac.lev)>1){
  for(i in 2:length(fac.lev)){
    df1=rbind(df1,
              data.frame(alp,
                         size,
                         ID=nrow(X)+i,
                         facet=fac.lev[i],
                         Group='load',
                         ppm=ppm[idx],
                         Intensity=cv1,
                         load=cols))
  }}

theme_bare <- theme(
  legend.position = "none",
  panel.background = element_blank(),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   plot.background = element_blank()
)


# create plots g=loading, g1=spec overlay
g=ggplot()+
  geom_line(data=df1, aes_string('ppm', 'Intensity', color='load', group='ID'), size=0.8)+
  scale_x_reverse(breaks=round(seq(shift[1], shift[2], by=abs(diff(shift))/20),3))+
  scale_y_continuous(limits=limY)+
  ggtitle(title)+
  xlab(expression(delta~{}^1*H~'(ppm)'))+
  ylab('Intensity (AU)')+
  facet_grid(facet~.)+
  theme_bare+
  theme(axis.text = element_text(colour="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #legend.position="none",
        legend.justification=c(1,1), legend.position=c(0.98,0.99))

if(type=='Statistical reconstruction'){
  g=g+scale_colour_gradientn(colors=matlab.like2(10), name='cor(t,x)', limits=raCol)
}else{
  g=g+scale_colour_gradientn(colors=matlab.like2(10), name=expression(abs~w[pred*','~sc]), limits=raCol)
}


g1=ggplot()+
  scale_x_reverse(breaks=round(seq(shift[1], shift[2], by=abs(diff(shift))/20),3))+
  scale_y_continuous(limits=limY)+
  ggtitle(title)+
  xlab(expression(delta~{}^1*H~'(ppm)'))+
  ylab('Intensity (AU)')+
  facet_grid(facet~.)+
  theme_bw()+
  theme(axis.text = element_text(colour="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(colour="white"),
        #legend.position="none",
        legend.justification=c(1,1), legend.position=c(0.98,0.55))+
  scale_colour_brewer(palette = "Dark2")


if(length(an)==1){
  g1=g1+
    geom_line(data=df, aes_string('variable', 'value', group='ID'), colour='black', alpha=alp, size=size)
}
if(length(an)==2){
  g1=g1+
    geom_line(data=df, aes_string('variable', 'value',group='ID', colour=names(an)[2]), alpha=alp, size=size)
}
if(length(an)==3){
  if(length(an[[2]])==1){
    g1=g1+
      geom_line(data=df, aes_string('variable', 'value',group='ID', linetype=names(an)[3]), alpha=alp, size=size, colour='black')
  }else{
    g1=g1+
      geom_line(data=df, aes_string('variable', 'value',group='ID', colour=names(an)[2] , linetype=names(an)[3]), alpha=alp, size=size)
  }
}




grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( g1 , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( g , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )

}
