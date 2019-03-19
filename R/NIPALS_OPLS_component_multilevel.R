#' Calculate PLS component, orthogonalise X and remove orthogonal information
#' @description Calculate PLS component, orthogonalise X and remove orthogonal information. This function is not for general use, rather than part of the opls function.
#' @param X Input matrix with rows and columns representing observations and variables
#' @param Y Dependend variable, in form of dummy matrix (multi-levels allowed) or numeric column vector
#' @return Returned is a list with the following entries:
#' \item{Filtered X}{Orthogonal filtered X matrix.}
#' \item{Scores X pred}{PLS component scores.}
#' \item{Loadings X pred}{PLS component loadings for X.}
#' \item{Weights pred}{PLS component variable weights for X.}
#' \item{Scores X orth}{Orthogonal component scores.}
#' \item{Loaadings X orth}{Orthogonal component loadings for X.}
#' \item{Weights X orth}{Orthogonal component X variable weights.}
#' \item{Loadings Y}{PLS component Y loadings.}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @noRd
NIPALS_OPLS_component_mulitlevel <- function(X, Y) {
    if (ncol(Y) > 1) {
        # 4 initialise scores u with column of Y, this is for two level outcome
        W_h <- apply(Y, 2, function(i, x = X) {
            (t(i) %*% x)/drop(crossprod(i))
        })
        ss_T <- ss_W_h <- sum(W_h^2)
        # estimate PCA the principal components of W as long as the ratio of SS of current score vector t divided by SS of W is larger than given
        # threshold, typicaly 10e-10
        count <- 1
        while ((ss_T/ss_W_h) > 1e-09) {
            # print((ss_T/ss_W_h))
            if (count == 1) {
                w_pca <- NIPALS_PCAcomponent(W_h)
                T_w <- rbind(w_pca$Scores)
            } else {
                w_pca <- NIPALS_PCAcomponent(w_pca$`Residual X`)
                T_w <- cbind(T_w, w_pca$Scores)
            }
            ss_T <- sum(w_pca$Scores[, 1]^2)
            count <- count + 1
            if (count > 30) {
                stop("endless loop?")
            }
        }
    }
    dd <- 1
    count <- 1
    u <- cbind(Y[, 1])
    while (dd > 1e-10) {
        # 1 calc weights scores and loadings
        w_h <- (t(u) %*% X)/drop(crossprod(u))
        # 2 normalise with norm
        w_h <- t(w_h)/sqrt(sum(w_h[1, ]^2))  # normalisation
        # 3 calc scores X
        t_h <- (X %*% w_h)/crossprod(w_h)[1, 1]
        # 4 calc loadings Y -> c is q
        c_h <- (t(t_h) %*% Y)/drop(crossprod(t_h))
        # 5. cal new u and compare with one in previus iteration (stop criterion)
        u_new <- (Y %*% t(c_h))/drop(crossprod(t(c_h)))
        if (count > 1) {
            dd <- sum((u_new - u)^2)/sum(as.numeric(u_new)^2)  #print(dd)
        }
        u <- u_new
        count <- count + 1
    }
    # 6 calc loadings
    p_h <- (t(t_h) %*% X)/drop(crossprod(t_h))
    # for multicolumn Y: estimate orthogonal component
    if (ncol(Y) > 1) {
        # 11: orthogonalise p_h (use the PCA scores of Y weights: T_w)
        for (i in 1:ncol(T_w)) {
            t_w <- cbind(T_w[, i])
            p_h <- p_h - ((t(t_w) %*% t(p_h)/crossprod(t_w)[1, 1]) %*% t(t_w))
        }
        w_o <- p_h
    } else {
        w_o <- p_h - ((t(w_h) %*% t(p_h)/drop(crossprod(w_h))) %*% t(w_h))
    }
    # 8: normalise
    w_o <- w_o/sqrt(sum(w_o[1, ]^2))  # normalisation
    # 9: orthogonal scores
    t_o <- (X %*% t(w_o))/drop(crossprod(t(w_o)))
    # 10 orthogonal laodings
    p_o <- (t(t_o) %*% X)/drop(crossprod(t_o))
    # 11: Filter data
    E_opls <- X - (t_o %*% p_o)
    # Y_res = Y - u %*% c_h 16: return paramters
    res <- list(`Filtered X` = E_opls, `Scores X pred` = t_h, `Loadings X pred` = p_h, `Weights pred` = w_h, `Scores X orth` = t_o, `Loadings X orth` = p_o, 
        `Weights X orth` = w_o, `Loadings Y` = c_h)
    return(res)
}
