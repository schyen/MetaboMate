#' An S4 class to represent OPLS models
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
# define slots for OPLS_Torben object
setClass("OPLS_MetaboMate", representation(type = "character", t_pred = "matrix", p_pred = "matrix", w_pred = "matrix", betas_pred = "numeric", 
    Qpc = "matrix", t_cv = "matrix", t_orth_cv = "matrix", t_orth = "matrix", p_orth = "matrix", w_orth = "matrix", nPC = "numeric", summary = "data.frame", 
    X_orth = "matrix", Y_res = "matrix", Xcenter = "numeric", Xscale = "numeric", Ycenter = "numeric", Yscale = "numeric", Yout = "data.frame", 
    Parameters = "data.frame", E = "matrix"))
