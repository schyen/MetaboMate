#' Predicting scores and Y outcome using PLS component(s)
#' @param X Input matrix with rows and columns representing new observations with identical number of variables than length of \code{pls_mod$loadings}.
#' @param pls_mod PLS component as outputted by function \cite{\code{NIPALS_PLS_component}}
#' @return Returned is a list with the following entries:
#' \item{Scores_pred}{Predicted PLS component scores.}
#' \item{Y_hat}{Predicted outcome. Averaged over all components in multi-PLS-component case.}
#' #' \item{Y_hat_components}{Predicted outcome (Y) for each component.}

  pls_prediction=function(pls_mod, X){

  # in case of one sample scenario
  if(is.null(ncol(X))){X=rbind(X)}

  # calc scores predictions and residuals
  t_pred = X %*% t(pls_mod$weights)
  E_h = X - (t_pred %*% pls_mod$loadings)

  betas=pls_mod$betas
  q_h=pls_mod$Qpc
  #if(q_h[1,1]==-q_h[2,1]){q_h=rbind(q_h[1,])}

  res=matrix(NA, nrow=nrow(X), ncol=ncol(pls_mod$scores))
  for(i in 1:ncol(pls_mod$scores)){
    opts=t(cbind(betas[i]) %*% t_pred[,i]) %*% (q_h [,i])
    res[,i]=apply(opts, 1, sum) # this sums up over all components, in two case scenario both columns are identical but different sign: that nulls this.
  }
  totalPrediction=apply(res, 1, sum)

  res=list('Scores_pred'=t_pred, 'Y_hat'=totalPrediction, 'Y_hat_components'=res)
  return(res)
}

