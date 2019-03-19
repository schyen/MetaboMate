#' Perform Principal Component Analysis
#' @export
#' @description This function is used to perform Principal Component Analysis (PCA).
#' @param X Numeric input matrix with each row representing an observation and each column a metabolic feature.
#' @param pc Desired number of principal components.
#' @param scale Desired scaling method: \code{None}, \code{UV} (unit variance) or \code{Pareto} (Pareto scaling).
#' @param method Algorithm for computing PCA. NIPALS is standard and usually fine. It can handle small amounts of missing/NA values.
#' @details Other methods include: 'svd', 'rnipals', 'bpca', 'ppca', 'svdImpute', 'robustPca', 'nlpca', 'llsImpute', 'llsImputeAll'. If these methods are specified, the \code{pca} function from the \code{pcaMethods} package is used to fit PCA model (see References).
#' @param center Logical indicating if data should be mean centered.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns a \emph{PCA_MetaboMate} S4 object.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @seealso \code{\link[=PCA_MetaboMate-class]{PCA_MetaboMate}} \code{\link[pcaMethods]{pca}} \code{\link{plotscores}} \code{\link{plotload}} \code{\link{opls}}
#' @importFrom pcaMethods pca
#'
pca <- function(X, pc = 2, scale = "UV", center = T, method = "nipals") {
    X <- as.matrix(center_scale(X, idc = "all", center, scale))
    if (method == "nipals") {
        res <- list()
        for (i in 1:pc) {
            if (i == 1) {
                res[[i]] <- NIPALS_PCAcomponent(X)
            } else (res[[i]] <- NIPALS_PCAcomponent(X = res[[i - 1]][[1]]))
        }
        Tpc <- sapply(res, "[[", 2)
        Ppc <- sapply(res, "[[", 3)
        # total.var<-sum(diag(cov(X))) #Calculate total variance in
        total.var <- totSS(X)
        prop.var <- rep(NA, ncol(Tpc))
        print(dim(Tpc))
        print(dim(Ppc))
        cum.var <- rep(NA, ncol(Tpc))  #Create #Calculate proportion of variance explained and cumulative
        for (i in 1:ncol(Tpc)) {
            prop.var[i] <- totSS(Tpc[, i] %o% Ppc[, i])/total.var
            # prop.var[i]<-var(Tpc[,i])/total.var
        }
        mod_pca <- new("PCA_MetaboMate", t = Tpc, p = Ppc, nc = pc, R2 = prop.var)
    } else {
        mod <- pcaMethods::pca(X, nPcs = pc, scale = "none", center = F, method = method)
        r2 <- mod@R2cum
        for (i in 1:pc) {
            if (i == 1) {
                next
            } else {
                r2[i] <- r2[i] - cumsum(r2[1:(i - 1)])
            }
        }
        mod_pca <- new("PCA_MetaboMate", t = mod@scores, p = mod@loadings, nc = pc, R2 = r2)
    }
    return(mod_pca)
}
