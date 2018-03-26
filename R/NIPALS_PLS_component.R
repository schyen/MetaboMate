#' Calculating a single PLS component
#' @param X Input matrix with rows and columns representing observations and variables
#' @param Y Dependend variable, in form of dummy matrix (multi-levels allowed) or numeric column vector
#' @return Returned is a list with the following entries:
#' \item{X.res}{Residual X matrix.}
#' \item{Y.res}{Residual Y matrix.}
#' \item{scores}{PLS component scores.}
#' \item{loadings}{PLS component loadings.}
#' \item{weights}{PLS component weights.}
#' \item{betas}{PLS X coefficients.}
#' \item{Q.pc}{PLS Y coefficient.}


NIPALS_PLS_component<-function(X, Y){

  # X=T*P +E
  # Y=U*Q +F*

  dd<-1
  count<-1

  # initialise scores (t_h)
  u_h <- cbind(Y[,1]) # u_h start is Y
  while(dd>1e-10){

    # X block: calc weights, normalise and cacultate scores (t_h)
    w_h <- t(X) %*% u_h / drop(crossprod(u_h))
    w_h <- w_h / sqrt(sum(w_h^2)) # scale to length 1

    # calc X scores (t_h)
    t_h <- X %*% w_h / drop(crossprod(w_h)) # normalisation not needed (is one)

    # check convergence of t with last iteration (dd crit in while loop)
    if(count>1){
      dd<-sum((t_h[,1]-t_old[,1])^2);
    }
    t_old<-t_h

    ## Y block, calc weights, normalise and calc Y scores
    ## steps can be omitted for 2 class Y (simply by setting q_h=1)
    q_h<-t(Y) %*% t_h /  drop(crossprod(t_h))
    q_h<-q_h/sqrt(sum(q_h^2)) # normalisation

    u_h<-(Y) %*% (q_h) / drop(crossprod((q_h)))
    count<-count+1
  }

    # calculate the X loadings and rescale scores and weights accordingly
    p_h <- (t(X) %*% t_h) / drop(crossprod(t_h))
    p_h <- p_h / sqrt(sum(q_h^2))

    # calc y loadings
    yl<-t(Y) %*% u_h/drop(crossprod(u_h))

    # calc b for residuals Y
    # calculate beta (regression coefficient) via inverse insted of subtracting q_h directly
    b<-(t(u_h) %*% t_h) /drop(crossprod(t_h))

    # # calc b for residuals Y
    # # calculate beta (regression coefficient) via inverse insted of subtracting q_h directly
    # b= solve(t(t_h) %*% t_h) %*% t(t_h) %*% cbind(u_h)
    # Y_res = Y - t( b[1,1] *t(t_h)) %*% (q_h)

  # calc residual matrice X and Y
  X_res <- X - (t_h %*% t(p_h))
  Y_res <- Y - (drop(b) *t_h) %*% t(q_h)


  # define slots for PLS_Torben object
  res<-list(X.res=X_res, Y.res=Y_res, scores=t_h, loadings=t(p_h), weights=t(w_h), betas=as.numeric(b), Qpc=q_h)
  #
  # track = setClass('PLS_Torben', slots=c(Xres='matrix', Yres='matrix', scores='matrix', loadings='matrix', weights='matrix', betas='numeric', Qpc='matrix')) #, Xcenter='numeric', Xscale='numeric', Ycenter='numeric', Yscale='numeric'
  #
  # # create slots for PLS_Torben objects
  # mod_pls=track(Xres=X_res, Yres=Y_res, scores=t_h, loadings=t(p_h), weights=t(w_h), betas=as.numeric(b), Qpc=q_h) #, Xcenter=meanX, Xscale=sdX, Ycenter=meanY, Yscale=sdY
  #

  return(res)
}

