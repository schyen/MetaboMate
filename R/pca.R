#' Perform Principal Component Analysis
#' @export
#' @description This function is used to perform Principal Component Analysis (PCA).
#' @param X Numeric input matrix with each row representing an observation and each column a metabolic feature.
#' @param pc Desired number of principal components.
#' @param scale Desired scaling, currently only unit variance scaling (or no scaling \code{scale=FALSE} is implemented.
#' @param center Logical indicating if data should be mean centred.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns a \emph{PCA_MetaboMate} S4 object.
#' @author Torben Kimhofer

pca=function(X, pc=2, scale='UV', center=T){

  X=as.matrix(MetaboMate:::center_scale(X, idc = 'all', center, scale))

  res=list()
  for(i in 1:pc){
    if(i==1){
      res[[i]]=MetaboMate:::NIPALS_PCAcomponent(X)
    }else(
      res[[i]]=MetaboMate:::NIPALS_PCAcomponent(X=res[[i-1]][[1]])
    )
  }

  Tpc=sapply(res, '[[',2)
  Ppc=sapply(res, '[[',3)

  #total.var<-sum(diag(cov(X))) #Calculate total variance in
  total.var<-totSS(X)
  prop.var<-rep(NA,ncol(Tpc));
  cum.var<-rep(NA,ncol(Tpc)) #Create #Calculate proportion of variance explained and cumulative
  for(i in 1:ncol(Tpc)){prop.var[i]<-var(Tpc[,i])/total.var}

  mod_pca=new('PCA_MetaboMate',
              t=Tpc,
              p=Ppc,
              nc=pc,
              R2=prop.var)

  return(mod_pca)

}




