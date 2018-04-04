#' Plotting overlayed NMR specra
#' @aliases specOverlay
#' @export
#' @param X Input NMR data matrix with row representing spectra.
#' @param ppm ppm vector with its length equals to \code{nrow(X)}.
#' @param shift Chemical shift area to be plotted. This should be kept as small as possible (see Details).
#' @param an List with one to three elements specifying facetting, colour and linetype (see Details).
#' @param alp Alpha value for lines (number between 0 and 1 whereas 0 is fully transparent and 1 is fully opaque).
#' @param title Plot title.
#' @param size Line width (0.5 is a good start).
#' @param ... Additional paramters passed on to ggplot's facet function.
#' @description  Plotting overlayed NMR specra. This function is based on ggplot2, a high-level plotting R package. For high ppm ranges computation time is relatively, so the range of input argument \code{shift} should be as small as possible. List argument \code{an} must have the first element define, even if it is only a single valuess. If colour and line width is specified, then at least one list elements of \code{an} must have the same length as \code{X}.
#' @seealso \code{\link{readBruker}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes_string scale_x_reverse ggtitle xlab facet_grid theme_bw theme element_text geom_line
#' @importFrom colorRamps matlab.like2
#' @importFrom scales pretty_breaks
#' @importFrom stats as.formula
#' @seealso \code{\link{facet_grid}}

specOverlay=function(X, ppm, shift=c(-0.01,0.01), an=list('facet' , 'col', 'ltype'), alp=0.7, size=0.5, title='', ...){

  if(is.null(names(an))){cat('No facet, colour and linetype names given. See an argument in ?specOverlay\n')
    names(an)=paste('an', 1:length(an), sep='')}

  le.arg<-paste(length(an))

  if ('' %in% names(an)){
    idx=which(names(an)=='')
    names(an)[idx]=paste('an', idx, sep='')
  }
  names(an)=gsub(' ', '.', names(an))


  idx=get.idx(shift, ppm)
  specs=X[,idx]
  colnames(specs)=paste("Idx", idx, sep='_')

  # if(le.arg==2){
  #   col.cat=is.factor(an[[2]])|is.character(an[[2]])|is.logical(an[[2]])
  #   }
  # if(le.arg==3){
  #   col.cat=is.factor(an[[2]])|is.character(an[[2]])|is.logical(an[[2]])
  #   an[[3]]=factor(an[[3]])
  #   }

  # create dataframe for ggplot function
  df=data.frame(do.call(cbind.data.frame, an), ID=1:nrow(specs), alp, specs)
  colnames(df)[1:le.arg]=names(an)
  df=melt(df, id.vars=c('alp','ID',names(an)))
  df$variable=ppm[as.numeric(gsub('Idx_', '', df$variable))]

  # initiate generic ggplot object
  g<- ggplot()+
    scale_x_reverse(breaks=seq(shift[1], shift[2], by=abs(diff(shift))/20),
                    name=expression(delta~{}^1*H~'(ppm)'))+
    scale_y_continuous(breaks=pretty_breaks(), name='Intensity (AU)')+
    ggtitle(title)+
    facet_grid(as.formula(paste(names(an)[1],"~ .")))+ #, ...
    theme_bw()+
    theme(axis.text = element_text(colour="black"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # add colour and line type
  switch(le.arg,
         '1' = {g<-g+
           geom_line(data=df, aes_string(x='variable',
                                         y='value',
                                         group='ID'),
                     colour='black',
                     alpha=alp,
                     size=size)
         },
         '2' = {g<-g+
               geom_line(data=df, aes_string(x='variable',
                                             y='value',
                                             group='ID',
                                             colour=names(an)[2]),
                         alpha=alp,
                         size=size)
         # add multi-colour gradient if colour vector is not factor/char
         col.cat=is.factor(an[[2]])|is.character(an[[2]])|is.logical(an[[2]])
         if(col.cat==F){
           g<-g + scale_colour_gradientn(colors=matlab.like2(length(an[[2]])))
          }
         },
         '3' = {an[[3]]=factor(an[[3]])
           g<- g+
               geom_line(data=df, aes_string(x='variable',
                                             y='value',
                                             group='ID',
                                             colour=names(an)[2],
                                             linetype=names(an)[3]),
                         alpha=alp, size=size)
         # add multi-colour gradient if colour vector is not factor/char
         col.cat=is.factor(an[[2]])|is.character(an[[2]])|is.logical(an[[2]])
         if(col.cat==F){
           g<-g + scale_colour_gradientn(colors=matlab.like2(length(an[[2]])))
          }
         }
         )


  return(g)
}
