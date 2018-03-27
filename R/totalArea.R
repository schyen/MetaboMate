#' Total area normalisation
#' @export
#' @aliases totalArea
#' @param X Input matrix where rows represent samples.
#' @param add.DilF Character string for new variable for dilution factor (will be exported as global envrionment variable). Can be left NULL if dilution factor should not be exported.
#' @author Torben Kimhofer
totalArea=function(X, add.DilF='ta.dilf'){
  X[X==0]=0.0001

  dil.F=1e10/apply(X, 1, sum)

  X.ta=t(apply(rbind(1:nrow(X)), 2, function(i){
    X[i,]*dil.F[i]
  }))

  if(!is.null(add.DilF)){
    assign(add.DilF, dil.F, envir=.GlobalEnv)
  }

  return(X.ta)
  }
