#' Generate k cross-validation sets
#' @param cv.k Number of cross validation sets (k-fold).
#' @param Y Dependent variable in vector or dummy matrix format (see Details).
#' @param type Indicating if a regression (\code{R}) or discriminant analysis (\code{DA}) will be performed.
#' @param method Type of cross validation: k-fold or Monte Carlo-Cross Validation (MCCV), sampling can be performed totally random or group stratified: 'k-fold', 'k-fold_stratified', 'MC', 'MC_stratified'.
#' @details If input argument \code{Y} is a dummy marix then function \cite{\code{centre_scale}} has been applied beforehand.The number of cross validation sets is related to the number of samples in each group. If in doubt, set \code{k=nrow(Y)} or \code{k=length(Y)} if Y is matrix or vector, respectively.
#' @seealso center_scale
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
cv_sets_method <- function(cv.k = 7, Y, type = c("R", "DA"), method = "k-fold_stratified") {
    if (type == "R" & (method == "k-fold_stratified" | method == "MC_stratified")) {
        cat("Stratified k-fold CV sampling not implemented for numeric Y. Setting opls function argument cv.type = \"k-fold\" and continue.\n")
        method <- "kfold"
    }
    if (is.null(ncol(Y))) {
        Y <- as.matrix(Y, ncol = 1)
    }
    switch(method, `k-fold` = {
        sets <- sample(1:cv.k, nrow(Y), replace = T)
        if (cv.k == nrow(Y)) {
            # n-fold case
            sets <- 1:nrow(Y)
        }
        sets_list <- lapply(1:cv.k, function(i, idc = sets) {
            which(idc != i)
        })
    }, MC = {
        sets_list <- lapply(1:cv.k, function(i, le = nrow(Y)) {
            sample(1:nrow(Y), nrow(Y) * 2/3)
        })
    }, `k-fold_stratified` = {
        levs <- as.vector(unique(Y))  # preserves the number of percentage of each class in each fold
        ct <- table(Y)
        nobs <- floor(min(ct)/cv.k)
        if (min(ct) <= cv.k) {
            cat("Number of observations in group ", names(which.min(ct)), " is too low (", min(ct), ")!\n", sep = "")
            return(NULL)
        }
        # cat('Minimum number of observations in each fold: ', names(which.min(ct)), ' = ', nobs, '\n', sep='')
        set_lev <- lapply(levs, function(x, k = cv.k, y = Y, nob = nobs) {
            idx <- sample(which(y == x))
            k1 <- rep(1:k, length.out = length(idx))
            names(k1) <- c(idx)
            return(k1)
        })
        sets_list <- lapply(1:cv.k, function(i, idc = unlist(set_lev)) {
            idx <- which(idc != i)
            as.numeric(names(idc)[idx])
        })
    }, MC_stratified = {
        levs <- as.vector(unique(Y))  # sample in a way that 2/3 training set, with equal number of samples per class in each fold
        ct <- table(Y)
        nobs <- floor(min(ct) * 2/3)
        sets_list <- lapply(1:cv.k, function(i, y = Y, nob = nobs, ctab = ct) {
            print(i)
            k1 <- list()
            for (j in 1:length(ctab)) {
                k1[[j]] <- sample(which(y == names(ctab)[j]), nob)
            }
            unlist(k1)
        })
    })
    return(sets_list)
}
