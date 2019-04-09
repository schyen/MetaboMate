#' Calculating a single OPLS component
#' @param X Input matrix with rows and columns representing observations and variables
#' @param Y Dependent variable, in form of dummy matrix (multi-levels allowed) or numeric column vector
#' @return Returned is a list with the following entries:
#' \item{Filtered X}{Orthogonal filtered X matrix.}
#' \item{Scores X pred}{PLS component scores.}
#' \item{Loadings X pred}{PLS component loadings for X.}
#' \item{Weights pred}{PLS component variable weights for X.}
#' \item{Scores X orth}{Orthogonal component scores.}
#' \item{Loaadings X orth}{Orthogonal component loadings for X.}
#' \item{Weights X orth}{Orthogonal component X variable weights.}
#' \item{Loadings Y}{PLS component Y loadings.}
#' @description Calculating a single OPLS component. This function should only be used for two level-Y, coded as a numeric vector.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @noRd

NIPALS_OPLS_component<-function(X, Y){

  # 4 initialise scores u with column of Y, this is for two level outcome
  u <- cbind(Y[,1])

  dd<-1
  count<-1
  while(dd>1e-10){
    # 1 calc weights scores and loadings
    w_h <- (t(u) %*% X )/ drop(crossprod(u))

    # 2 normalise with norm
    w_h <- t(w_h) / sqrt(sum(w_h[1,]^2)) # normalisation

    # 3 calc scores X
    t_h <- (X %*% w_h) / crossprod(w_h)[1,1]

    # if(count>1){
    #   dd=sum((t_h[,1]-t_old[,1])^2); #print(dd)
    # }

    # 4 calc loadings Y -> c is q
    c_h<-(t(t_h) %*% Y) /  drop(crossprod(t_h))

    #c_h=c_h/sqrt(sum(c_h[1,]^2)) # normalisation
    # 5. cal new u and compare with one in previus iteration (stop criterion)
    u_new<-(Y %*% t(c_h)) / drop(crossprod(t(c_h)))

    if(count>1){
      dd<-sum((u_new-u)^2)/sum(as.numeric(u_new)^2); print(dd)
    }
    u<-u_new
    count<-count+1
  }
  # 6 calc loadings
  p_h <- (t(t_h) %*% X) / drop(crossprod(t_h))


  # # for multicolumn Y: estimate orthogonal component
  # # 11: orthogonalise p_h (use the scores of Y weights: T_w)
  # for(i in 1:ncol(T_w)){
  #   t_w=cbind(T_w[,i])
  #   p_h=p_h - ((t(t_w) %*% t(p_h)/crossprod(t_w)[1,1]) %*% t(t_w))
  # }
  #
  #
  # step 7
  w_o<-p_h - ((t(w_h) %*% t(p_h)/drop(crossprod(w_h))) %*% t(w_h))

  # 8: normalise
  w_o<-w_o/sqrt(sum(w_o[1,]^2)) # normalisation

  # 9: orthogonal scores
  t_o<-(X%*% t(w_o)) / drop(crossprod(t(w_o)))

  # 10 orthogonal laodings
  p_o<-(t(t_o) %*% X) / drop(crossprod(t_o))



  # 16: return paramters
  # 11: Filter data
  E_opls <- X - (t_o %*% p_o)
  #Y_res = Y - u %*% c_h
  res<-list('Filtered X'=E_opls, 'Scores X pred'=t_h, 'Loadings X pred'=p_h, 'Weights pred'=w_h, 'Scores X orth'=t_o, 'Loadings X orth'=p_o, 'Weights X orth'=w_o,'Loadings Y'=c_h)

  return(res)
}




