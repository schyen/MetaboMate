# NIPALS PCA algorithm
# Partial least squares and regression: a tutorial (Geladi and Kowalski, Analytical Chimica Acta, 1986)

#' Calculating a single PCA component
#' @param X Input matrix with rows and columns representing observations and variables.
#' @return Returned is a list with the following entries:
#' \item{Residual X}{Residual X matrix.}
#' \item{Scores}{PCA component scores.}
#' \item{Loadings}{PLS component loadings.}
#' @references {Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. Analytica Chimica Acta 185, 1-17.}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
NIPALS_PCAcomponent<-function(X){
  # initialise scores (t_h)
  t_h<-cbind(X[,1])
  count<-1
  dd<-5
  while(dd>1e-5){
    # calc first loadings (p_h)
    p_h<- t(t_h) %*% X / crossprod(t_h)[1,1]

    #normalise p_h to length 1
    p_h<- p_h /sqrt(sum(p_h[1,]^2))

    # update scores and repeat until difference is small
    t_new<-(X %*% t(p_h)) / crossprod(t(p_h))[1,1]

    dd<-sum((t_h[1,]-t_new[1,])^2); #print(dd)

    t_h<-t_new
    if(count>5000){stop('Failed to converge!!!')}
    #print(count)
    count<-count+1
  }
  cat('iterations: ', count, '\n', sep='')

  X_res<-X-(t_h %*% (p_h))

  res<-list('Residual X'=X_res, 'Scores'=t_new, 'Loadings'=p_h)
  return(res)
}
