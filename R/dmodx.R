#' Calculating distance to the model in X space
#' @export
#' @param model OPLS model of type \code{OPLS_MetaboMate}.
#' @param plot Logical indicating if results should be visualised.
#' @return The projection distance of each observation in the model (\code{DModX}).
#' @references Bylesjö M., \emph{et al.} (2002) OPLS discriminant analysis: combining the strengthsof PLS-DA and SIMCA classification. \emph{Journal of Chemometrics}, 20, 341-51.
#' @references Wold S. (1976) Pattern recognition by means of disjoint princi-pal components models.  \emph{Pattern Recognition}, 8, 127-39.
#' @importFrom ggplot2 ggplot aes_string geom_point scale_colour_gradientn geom_hline xlab scale_y_continuous theme_bw theme element_blank element_text
#' @importFrom colorRamps matlab.like
#' @importFrom stats t.test sd
# E=residual Matrix
# N=number of samples
# K=number of variables
# A=number of model components
# A0= (1 if mean centred, 0 otherwise)

dmodx=function(model, plot=T){

  if(class(model)!='OPLS_MetaboMate'){stop('Please provide a OPLS_MetaboMate object.')}

  E=model@E
  N=nrow(E)
  K=ncol(E)
  A=ncol(model@t_pred) # in case of OPLS-DA (alwasy one predictive component)

  if(as.logical(model@Parameters$Value[model@Parameters$Paramter=='Center'])==T){A0=1}else{A0=0}
  # loop over all observations in residual matrix, calc SS residuals / observations and normalise by TSS
  ss_res=apply(E, 1, function(x) sum(x^2))
  dmodX=sqrt(ss_res/(K-A)) / sqrt(sum(ss_res)/((N-A-A0)*(K-A)))

  tt=t.test(dmodX, alternative='less')
  ci95=tt$conf.int[2]+2*sd(dmodX)
  df=data.frame(col=model@t_cv, ID=1:length(dmodX), DmodX=dmodX, passedT.test=dmodX<tt$conf.int[2]+2*sd(dmodX))
  if(plot==T){
    g=ggplot(df, aes_string('ID', 'DmodX', colour='col'))+
      geom_point()+
      scale_colour_gradientn(colours=matlab.like2(10), name=expression(t[pred]))+
      scale_y_continuous(limits = c(min(dmodX), max(c(dmodX, ci95))), name='DModX')+
      geom_hline(yintercept=ci95, linetype=2, colour='gray80')+
      xlab('Sample index')+
      theme_bw()+
      theme(panel.grid.minor.x = element_blank(),
            axis.text=element_text(colour="black"))
    plot(g)
  }

  return(df[,-1])
}


