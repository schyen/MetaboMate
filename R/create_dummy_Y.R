#' Create dummy matrix and check numeric variable for OPLS analysis
# @export
#' @param Y Dependent variable (numeric or factor for regression or discriminat analysis, resepctively)
#' @return List of two: Dummy matrix and data frame mapping Y-levels to numeric representations.
#' @aliases create_dummy_Y
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @noRd
create_dummy_Y <- function(Y) {
    if (!is.numeric(Y)) {
        Y_levels <- unique(Y)
        if (length(Y_levels) == 2) {
            Y_new <- cbind(as.numeric(as.factor(Y)))
        } else {
            Y_new <- matrix(-1, nrow = length(Y), ncol = length(Y_levels))
            for (i in 1:length(Y_levels)) {
                Y_new[which(Y == Y_levels[i]), i] <- 1
            }
            colnames(Y_new) <- Y_levels
        }
        Y_levs <- unique(data.frame(Original = Y, Numeric = Y_new, stringsAsFactors = F))
        # return(list(Y_new, Y_levs))
        return(list(Y_new, Y_levs))
    } else {
        return(list(cbind(Y), data.frame()))
    }
}
