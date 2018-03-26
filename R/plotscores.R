#' Plotting PCA, PLS or OPLS model scores
#' @export
#' @aliases plotscores
#' @description Function for plotting PCA, PLS or OPLS model scores.
#' @param model PCA, PLS or OPLS model of the package \emph{MetaboMate}.
#' @param pc Specifies which model components should be plotted on abcissa and ordinate (see Details).
#' @param an Colour and label specification (see Details).
#' @param title Plot title.
#' @param qc Row indices of QC samples. Can be left NA, however, if specified QC samples will be highlighted in the plot.
#' @param legend Position of the plot legend, set to NA if legend should be outside of the plotting area.
#' @param cv.scores Logical value (TRUE or FALSE) indicating if cross-validated scores should be plotted for the predictive component(s) (only for PLS and O-PLS).
#' @param ... Additional values passed on to \code{\link[ggplot2]{scale_colour_gradientn}}.
#' @details Scores colouring is specified with the argument \code{an}, which is a list of three elements. The first list element specifies the coulour (class factor required for categorical variables), the second list element specifies a point labeling (class character or factor) and the third list element specifies point shape. The Hotelling's \out{T<sup>2</sup>} ellipse is automatically included and calculated for the dimensions selected by the \code{pc} argument.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @references Hotelling, H. (1931) The generalization of Student’s ratio. \emph{Ann. Math. Stat.}, 2, 360-378.
#' @return This function returns a \emph{ggplot2} plot object.
#' @seealso \code{\link{MetaboMate-class}}
#' @seealso \code{\link{opls}}
#' @importFrom stats cov
#' @importFrom ellipse ellipse
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_polygon geom_hline geom_vline geom_point scale_colour_manual ggtitle labs theme labs scale_colour_gradientn scale_x_continuous
#' @importFrom colorRamps matlab.like
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics plot
#' @importFrom reshape2 melt
#' @importFrom scales pretty_breaks


plotscores=function(model, pc=c(1,2), an, title='', qc=NA, legend='in', cv.scores=T, ...){

  # define plot dimensions
  switch(class(model),
         'pcaRes' = {
           x <- model@"scores"[,pc[1]]
           y <- model@"scores"[,pc[2]]
         },
         "PCA_MetaboMate" = {
             x <- model@t[,pc[1]]
             y <- model@t[,pc[2]]
           },
         "OPLS_MetaboMate" = {
           if(cv.scores==T){
             x <- model@t_cv[,1]
             y <- model@t_orth_cv[,pc[1]]
           }else{x <- model@t_pred[,1]
           y <- model@t_orth[,pc[1]]
           }
         })

  ###### DATA PREP
  # prep df for ggplot2
  sc=data.frame(x,y, cl=an[[1]])
  if(length(an)==2){
    sc$shape=factor(an[[2]])}else{
      sc$shape=16}
  if(length(an)==3){
    sc$labs=an[[3]]}

  # prep sub-df for QC samples
  if(!is.na(qc[1])){
    df.qc=data.frame(x=x[qc],y=y[qc], cl=an[[1]][qc])
  }

  # Calculate Hotellings T2 elipse
  SD=cov(cbind(x,y))
  el=ellipse(SD, centre=colMeans(cbind(x,y)), level=0.95)
  colnames(el)=c("V1", "V2")

  # define axis ranges
  xlim=c(min((c(el[,1], x))), max((c(el[,1], x))))
  xlim=xlim+c(diff(range(xlim))*-0.05, diff(range(xlim))*+0.05)
  ylim=c(min((c(el[,2], y))), max((c(el[,2], y))))
  ylim=ylim+c(diff(range(ylim))*-0.05, diff(range(ylim))*+0.05)

  # prep df for Hotelltings T2 ellipse and axis label positions
  df=as.data.frame(el)
  y.labb=diff(range(ylim))*0.015
  x.labb=diff(range(xlim))*0.015
  df1=sc
  df1[,1]=df1[,1]+x.labb
  df1[,2]=df1[,2]+y.labb

  ###### COLOUR PREP
  # populate colours if number of levels exceeds colours in palette (n=8), fix point shape
  # minimum nb of factors for categorical variables in 17
  if(class(an[[1]])!='numeric' & length(table(an[[1]]))<17){
    type='categorical'
    if(length(unique(sc$cl))>8){
      cols=c(brewer.pal(8,"Dark2"), brewer.pal(8,"Set1")[1:(length(unique(sc$cl))-8)])}else{
        cols=brewer.pal(8,"Dark2")[1:length(unique(sc$cl))]}
  }else{
    type='continuous'
  }

  ###### FURTHER GGPLOT COMMANDS
  # prep ggplot2 skeleton
  g=ggplot()+
    geom_polygon(data=df, aes_string(x='V1', y='V2'), fill="white", alpha=0.4)+
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +
    ggtitle(title)+
    labs(colour=names(an)[1])

  # define colour schemes and create ggplot2 objects
  # if categorical colour vector has too many levels (>=17) groups will be converted to numeric and plotted with a colour scale
  # following if expression is for coloured points, no shape or label specifications
  if(length(an)==1){
   switch(type,
          'categorical'={
            g=g+
              geom_point(data=sc, aes_string('x','y', colour='cl'), alpha=0.7, shape=16)+
              scale_colour_manual(values=cols)+
              labs(colour=names(an)[1])
          },
          'continuous'={
            g=g+
              geom_point(data=sc, aes_string('x','y', colour='cl'), shape=16,  alpha=1)+#size=2.5,
              scale_colour_gradientn(colours=matlab.like(10),  breaks=pretty_breaks(), ...)+
              ggtitle(title)
          })
   }

  if(length(an)==2){
    switch(type,
           'categorical'={
             g=g+
               geom_point(data=sc, aes_string('x','y', colour='cl', shape='shape'), alpha=0.7)+
               scale_colour_manual(values=cols)+
               labs(colour=names(an)[1], shape=names(an)[2])
           },
           'continuous'={
             g=g+
               geom_point(data=sc, aes_string('x','y', colour='cl', shape='shape'), alpha=0.7)+
               scale_colour_gradientn(colours=matlab.like(10),  breaks=pretty_breaks(), ...)+
               ggtitle(title)
           })
    }

  if(length(an)==3){
    switch(type,
           'categorical'={
             if(length(unique(sc$shape))==1){
               g=g+
                 geom_point(data=sc, aes_string('x','y', colour='cl'), shape='shape', alpha=0.7)+
                 scale_colour_manual(values=cols)+
                 geom_text_repel(data=sc, aes_string('x','y', label='labs'))+
                 labs(colour=names(an)[1], shape=names(an)[2])
             }else{
               g=g+
                 geom_point(data=sc, aes_string('x','y', colour='cl', shape='shape'), alpha=0.7)+
                 scale_colour_manual(values=cols)+
                 geom_text_repel(data=sc, aes_string('x','y', label='labs'))+
                 labs(colour=names(an)[1], shape=names(an)[2])
             }
           },
           'continuous'={
             stop('A continuous variable can not be mapped to shape')
           })
  }

  if(!is.na(qc[1])){g=g+
    geom_point(data=df.qc, aes_string('x','y'), alpha=0.7, colour='black')}

  ###### AXIS LABELLING
  switch(class(model),
         'pcaRes' = {
           g=g+
             scale_x_continuous(name=paste("PC ", pc[1], " (",round(model@"R2"[pc[1]]*100,1)," %)", sep="" ))+
             scale_y_continuous(name=paste("PC ", pc[2], " (",round(model@"R2"[pc[2]]*100,1)," %)", sep="" ))
         },
         'PCA_MetaboMate' = {
           g=g+
             scale_x_continuous(name=paste("PC ", pc[1], " (",round(model@"R2"[pc[1]]*100,1)," %)", sep="" ))+
             scale_y_continuous(name=paste("PC ", pc[2], " (",round(model@"R2"[pc[2]]*100,1)," %)", sep="" ))
         },
         'OPLS_MetaboMate'= {
           if(ncol(model@t_orth)>1){
             comp='orthogonal components'}else{
               comp='orthogonal component'}
           if(cv.scores==F){
             if(model@type=='DA'){
               g=g+scale_x_continuous(name=expression(t[pred]))+
                 scale_y_continuous(name=expression(t[orth]))+
                 ggtitle(paste('OPLS - ', model@nPC," ",
                               comp,
                               " (R2X=",round(model@'summary'$'R2X'[model@nPC],2),
                               ', R2Y=',round(model@'summary'$'R2Y'[model@nPC],2),
                               ', Q2=',round(model@'summary'$'Q2'[model@nPC],2),
                               ', AUROC=',round(model@'summary'$'AUROC'[model@nPC],2),')'
                               , sep="" ))+
                 theme(plot.title = element_text(size=10))} else{
                 g=g+scale_x_continuous(name=expression(t[pred]))+
                   scale_y_continuous(name=expression(t[orth]))+
                   ggtitle(paste('OPLS - ', model@nPC," ",
                                 comp,
                                 " (R2X=",round(model@'summary'$'R2X'[model@nPC],2),
                                 ', R2Y=',round(model@'summary'$'R2Y'[model@nPC],2),
                                 ', Q2=',round(model@'summary'$'Q2'[model@nPC],2),')'
                                 , sep="" ))+
                   theme(plot.title = element_text(size=10))
               }
           }else{
             if(model@type=='DA'){
               g=g+
                 scale_x_continuous(name=expression(t[pred][','][cv]))+
                 scale_y_continuous(name=expression(t[orth][','][cv]))+
                 ggtitle(paste('OPLS - ',
                               model@nPC," ",
                               comp,
                               " (R2X=",round(model@'summary'$'R2X'[model@nPC],2),
                               ', R2Y=',round(model@'summary'$'R2Y'[model@nPC],2),
                               ', Q2=',round(model@'summary'$'Q2'[model@nPC],2),
                               ', AUROC=',round(model@'summary'$'AUROC'[model@nPC],2),')'
                               , sep="" ))+
                 theme(plot.title = element_text(size=10))
             }else{
               g=g+
                 scale_x_continuous(name=expression(t[pred][','][cv]))+
                 scale_y_continuous(name=expression(t[orth][','][cv]))+
                 ggtitle(paste('OPLS - ',
                               model@nPC," ",
                               comp,
                               " (R2X=",round(model@'summary'$'R2X'[model@nPC],2),
                               ', R2Y=',round(model@'summary'$'R2Y'[model@nPC],2),
                               ', Q2=',round(model@'summary'$'Q2'[model@nPC],2),')'
                               , sep="" ))+
                 theme(plot.title = element_text(size=10))
             }

           }
         })

  if(legend=='in'){
        g=g+
      theme(legend.position=c(1.01,1.02),legend.justification=c(1,1))}

  return(g)
}
