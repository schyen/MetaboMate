#' Calculating Total Sum of Squares (TSS)
#' @param X Mean centered and unit variance scaled matrix (see Details).
#' @return Total sum of squares.
#' @details Input matrix rows correspond to n observations and columns to m variables. For large matrices (here defined by m>3k), X is subsampled to speed up computation time. Currently, subsampling is achieved by selecting each selecting each third row element.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}

totSS<-function(X){
  # to speed calculations up, sample X matrix if large enought
  if(ncol(X) > 3000){ #nrow(X)<1000 |
    idx<-seq(from=1, to=ncol(X), by=3)
    tss<-sum(apply(X[,idx], 1, function(x){crossprod(x)}))
  }else{
    tss<-sum(apply(X, 1, function(x){crossprod(x)}))
  }
  return(tss)
}

