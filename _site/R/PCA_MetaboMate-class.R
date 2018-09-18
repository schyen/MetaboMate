#' An S4 class to represent PCA models
# define slots for PCA object
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
setClass('PCA_MetaboMate', representation(
    algorithm = 'character',
    t = 'matrix',
    p = 'matrix',
    nc = 'numeric',
    R2 = 'numeric'
  )
)

