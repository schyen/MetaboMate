#' Fitting Orthogonal-Partial Least Squares Models - parallelised
#' @export
#' @description This function is used to fit  Orthogonal-Partial Least Squares (O-PLS) models. It can be used to carry out regression or discriminant analysis. In the latter case the outcome can have two or more levels.
#' @param X Numeric input matrix (measurements derived by NMR spectroscopy or MS) with each row representing an observation and each column a metabolic feature.
#' @param Y Response vector or matrix with same length or number of columns than rows in X, respectively.
#' @param t_pred Parameter specifying the maximum number of predictive components (needed only for multi-factor Y)
#' @param center Logical value (TRUE or FALSE) indicating if features should be mean centered.
#' @param scale Desired scaling method (currently only no or unit variance scaling (UV) implemented).
#' @param cv.k The number of cross-validation sets. This depends on the number of observations in X but typically takes a value between 3 and 9.
#' @param cv.type Type or cross-validation: 'k-fold', 'k-fold_stratified', 'MC', 'MC_stratified' (see Details).
#' @param plotting Logical value (TRUE or FALSE) indicating if model parameters (R2X, Q2, etc) should be visualised once the model is trained.
#' @param maxPCo The maximum number of orthogonal components (in case stop criteria fail).
#' @details Models are fully statistically validated, currently only k-fold cross validation (CV) and class-balanced k-fold cross validation is implemented. Further extensions, e.g. Monte-Carlo CV, are work in progress. Although the algorithm accepts three and more levels as Y, model interpretation is more straightforward for pairwise group comparisons.
#' @references Trygg J. and Wold, S. (2002) Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16.3, 119-128.
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return This function returns an \emph{OPLS_MetaboMate} S4 object.
#' @seealso \code{\link{OPLS_MetaboMate-class}} \code{\link{dmodx}} \code{\link{plotscores}} \code{\link{plotload}} \code{\link{specload}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom graphics plot
#' @importFrom methods getSlots new representation setClass
#' @importFrom stats cov sd var
#' @importFrom utils methods
#' @importFrom ggplot2 ggplot aes aes_string scale_fill_manual scale_y_continuous theme_bw labs scale_x_discrete scale_alpha theme element_blank element_line element_rect geom_bar element_text
#' @importFrom pROC roc multiclass.roc
#' @importFrom reshape2 melt
#' @importFrom scales pretty_breaks
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom foreach foreach

opls.par <- function(X,
                                     Y,
                                     t_pred = 1,
                                     center = T,
                                     scale = 'UV',
                                     cv.k = 7,
                                     cv.type = 'k-fold',
                                     plotting = T,
                                     maxPCo = 5) {
  {

    cat('Setting-up parallel compouting ... ', sep = '')
    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)
    cat(no_cores, 'cores ...', sep = '')
    cat('done.\n', sep = '')

    # Determine if regression (R) or discriminant analysis (DA),
    # check Y for NA or non-numeric values
    Y1 <- Y
    if (is.numeric(Y)) {
      type <- 'R'
      if (is.infinite(max(abs(Y))) | anyNA(Y)) {
        stop('Y contains non-numeric values.')
      }
      Y_out <- create_dummy_Y(Y)
      Y <- Y_out[[1]]
    } else{
      type <- 'DA'
      if (anyNA(Y)) {
        stop('Y vector contains N/A\'s')
      }
      Y_out <- create_dummy_Y(Y)
      Y <- Y_out[[1]]
    }

    # check dimensions and convert to matrix, stop if features are non-numeric
    if (nrow(X) != nrow(Y)) {
      stop('X and Y dimensions do not match.')
    }
    X <- as.matrix(X)
    if (is.infinite(max(abs(X))) |
        anyNA(X)) {
      stop('X matrix contains non-numeric values.')
    }

    # Exclude features with sd=0
    # foreach(1:nrow(X))
    sdX <- apply(X, 2, sd)
    idx <- which(sdX == 0)
    if (length(idx) > 0) {
      cat(
        'A total of ',
        length(idx),
        ' features have a standard deviation of zero and will be excluded from analysis.\n',
        sep = ''
      )
      X <- X[, -idx]
      sdX <- sdX[-idx]
    }

    # Dats scaling and centering
    meanX <- apply(X, 2, mean)
    XcsTot <- center_scale(X,
                            idc = 'all',
                            center = meanX,
                            scale = sdX)
    YcsTot <- center_scale(Y,
                            idc = 'all',
                            center = T,
                            scale = 'UV')

    # Calculate Total Sum of Squares (TSS), required for calculation of Q2 etc
    # cat('Calculate total sum of squares...')
    tssx <- MetaboMate:::totSS(XcsTot)
    tssy <- totSS(YcsTot) # total variance that a model can explain (theoretical reference model, used for normalisation)
     cat('done.\n')
  }

  cat('Performing OPLS-', type, ' ... ', sep='')
  # k-fold cross-validation for calculation of Q2/AUROC
  # generate indices defining training set in each CV round
  cv_sets <- cv_sets_method(cv.k, Y=Y1, method = cv.type, type)

  # initialisations
  R2Y <- rss <- Press <- aucs <- rssx <- array()
  cv.res <- output <- list()
  t_cv <-  t_orth_cv <- preds <- matrix(NA, ncol = 1, nrow = nrow(X))
  nc <- 1
  enough <- F

  # cat('Fiting components...')

  # fit as many orthogoanl components until a stop criterion is TRUE
  while (enough == F) {

    # cat('Component ', nc, '...')
    if (nc>1) { # add column for nc > 1
       Xo=X
       Yo=Y
  }
        # fit each cv training set fit component and predict
        # validation set. Prediction accuracy output of is be used in stop criterion
    res=lapply(cv_sets, function(idc, Xcv=Xo, Ycv=Yo, ce=center, scale=sc){
      # scale/centre using CV trainine set samples
      idc <- cv_sets[[k]]
      Xcs <- center_scale(Xcv, idc, center=cen, scale=sc)
      Ycs <- center_scale(Ycv, idc, center=cen, scale=sc)

      # filter data: calc PLS component, then orthogonalise X with t_pred (=> t_o, p_o),
      E_opls <- NIPALS_OPLS_component_mulitlevel(X = Xcs[idc, ], Y = cbind(Ycs[idc, ]))

      p_orth=E_opls$`Loadings X orth`
      w_orth=E_opls$`Weights X orth`
      Xcv_res=E_opls$`Filtered X`

      # calculate predictive component with filtered matrix (Xres)
      pls_comp = NIPALS_PLS_component(X = cXcv_res, Y = cbind(Ycs[idc, ]))

      E_new_orth = Xcs[-idc, ]
      # Prediction of test data
      t_orth = E_new_orth %*% t(t(w_orth)) / drop(crossprod(t(t(w_orth))))
      e_new_orth = E_new_orth - (t_orth %*% t(p_orth))

      pred = pls_prediction(pls_mod = pls_comp, X = e_new_orth)

      return(list(preds, t_cv, t_orth_cv, ))
    })

    idc=order(unlist(cv_sets))
    preds=res[[1]][idc]
    t_cv=res[[2]][idc]
    t_orth_cv=res[[3]][idc]

    print(preds)
    stopImplicitCluster()

    # calculate rss, press and Q2
    Press[nc] = sum(apply(YcsTot, 2, function(x) {
      sum((x - preds[, nc]) ^ 2)
    })) / ncol(YcsTot)
    #print(Press)
    Q2_1 = 1 - (Press / tssy)
    #print(Q2_1)

    if (type == 'DA') {
      if (ncol(YcsTot) == 2) {
        mod = roc(response = Y[, 1], predictor = preds[, nc])
        aucs[nc] = mod$auc
      } else{
        mod = multiclass.roc(response = Y1, predictor = preds[, nc])
        aucs[nc] = mod$auc
      }


    } else{
      aucs[nc] = 0
    }


    # calculate r2Y with model using all data (no cv)
    if (nc == 1) {
      # initialisation
      # OPLS-filter X matrix repetitively until some stop criterion, save orthogonal loadings, scores and weights
      # need to scale this again

      E_opls_all = NIPALS_OPLS_component(X = XcsTot, Y = YcsTot)

      output[['t_orth']] = E_opls_all$`Scores X orth` # don't need scores for any calculation (only cv scores, and these are further down in prediction bit)
      output[['p_orth']] = E_opls_all$`Loadings X orth`
      output[['w_orth']] = E_opls_all$`Weights X orth`
      Xres = E_opls_all$`Filtered X`
    } else{
      E_opls_all = NIPALS_OPLS_component_mulitlevel(Xres, YcsTot)

      output[['t_orth']] = cbind(output[['t_orth']], E_opls_all$`Scores X orth`)
      output[['p_orth']] = rbind(output[['p_orth']], E_opls_all$`Loadings X orth`)
      output[['w_orth']] = rbind(output[['w_orth']], E_opls_all$`Weights X orth`)
      Xres = E_opls_all$`Filtered X`
    }

    # make predictions and calc R2
    # calculate predictive component with filtered data
    pls_comp_all = NIPALS_PLS_component(Xres, YcsTot)
    pred_all = pls_prediction(pls_mod = pls_comp_all, X = E_opls_all$`Filtered X`)
    #rss = sum((pred_all$Y_hat - YcsTot[, 1]) ^ 2)
    rss = sum(apply(YcsTot, 2, function(x) {
      sum((x - pred_all$Y_hat) ^ 2)
    })) / ncol(YcsTot)
    R2Y[nc] = 1 - (rss / tssy)

    X_ex = pls_comp_all$scores %*% pls_comp_all$loadings
    rssx[nc] = totSS(X_ex)
    R2x = (rssx / tssx)


    ### define stop criteria for fitting components
    if (length(which(is.na(Q2_1))) > 0) {
      cat('Something went wrong, Q2 is NA for component', nc, '\n',  sep = '')
    }

    if (nc == 1 & Q2_1[nc] < 0.05) {
      cat('At first PC, Q2 < 0.03: ', round(Q2_1[nc], 3), '\n',  sep = '')
      print('No sign. orthogonal components!')
      return(NULL)
      #print(Q2_cv)
      enough = T
    }


    if (nc > 1) {
      if ((Q2_1[nc] - Q2_1[nc - 1]) < 0.05) {
        # if q2 does not rise by more than 0.03 with new component
        # cat('At PC ', nc, ': delta Q2 < 0.03: ', round((Q2_1[nc] - Q2_1[nc -                                                                  1]), 3), '\n',  sep = '')
        enough = T
      }

    }

    if (Q2_1[nc] > 0.98) {
      # cat('At PC ', nc, ', Q2 > 0.98: ', round(Q2_1[nc], 3), '\n', sep = '')
      enough = T
    }

    if (nc == maxPCo) {
      enough = T
      nc = nc + 1
    }

    if (enough == F) {
      nc = nc + 1
    }

  }
  cat('done.\nA model with 1 predictive and', nc-1, 'orthogonal component(s) was fitted.\n\n')
  if (type == 'DA') {
    model.summary = data.frame(
      R2X = round(R2x, 2),
      R2Y = round(R2Y, 2),
      Q2 = round(Q2_1, 2),
      AUROC = round(aucs, 2)
    )
  } else{
    model.summary = data.frame(
      R2X = round(R2x, 2),
      R2Y = round(R2Y, 2),
      Q2 = round(Q2_1, 2)
    )
  }
  rownames(model.summary) = paste('PC_o', 1:(nrow(model.summary)))

  mm = cbind('PC_o' = rownames(model.summary), model.summary)
  mm = melt(mm, id.vars = 'PC_o')
  mm$PC = factor(mm$'PC_o', levels = rownames(model.summary))

  mm$alpha1 = 1
  mm$alpha1[mm$'PC_o' == paste('PC_o', nrow(model.summary))] = 0.7

  g = ggplot(mm, aes_string('PC_o',' value', fill = 'variable')) +
    geom_bar(stat = 'identity',
             position = "dodge",
             colour = NA,
             aes(alpha= mm$'alpha1') )+
    scale_y_continuous(limits = c(min(c(0, min(Q2_1) - 0.02)), 1), breaks = pretty_breaks()) +
    theme_bw() +
    labs(y = '',
      #title = paste('An  O-PLS-', type, ' model (1+', nc - 1, ') was fitted', sep = '')
      title = paste('O-PLS-', type, '  - Model Summary', sep = '')
    ) +
    scale_x_discrete(labels = paste('1+', 1:nc, sep = ''), name = 'Predictive + Orthogonal Component') +
    scale_fill_manual(
      values = c(
        'R2X' = 'lightgreen',
        'R2Y' = 'lightblue',
        'Q2' = 'red',
        'AUROC' = 'black'
      ),
      labels = c(
        expression(R ^ 2 * X),
        expression(R ^ {
          2
        } * Y),
        expression(Q ^ 2),
        expression(AUROC[cv])
      ),
      name = ''
    ) +
    scale_alpha(guide = F, limits=c(0,1)) +
    theme(
      legend.text.align = 0,
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = 'black', size =
                                          0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = 'black', size = 0.55),
      axis.line.y = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(colour = 'white'),
      text=element_text(family="Helvetica")
    )

  model.summary=model.summary[-nrow(model.summary),]
  #
  if (plotting == T) {
    plot(g)
  }

  if (nc > 0) {
    # re calc t_pred with number of nc-2 (one componenet remove and nc counter has already been set 1 up)
    Xorth = cbind(output[['t_orth']][, 1:(nc - 1)]) %*% rbind(output[['p_orth']][1:(nc -
                                                                                      1), ])
    Xres = XcsTot - Xorth
    pls_comp_all = NIPALS_PLS_component(Xres, YcsTot)

    E = Xres - (pls_comp_all$scores %*% pls_comp_all$loadings) # this is for calculation of DModX

    # this is for visualisation of orthogonal variation (currently not implemented)
    #E_pca=NIPALS_PCAcomponent(E+Xorth)



    dd = data.frame(
      Paramter = c(
        'Center',
        'Scale',
        'nPred',
        'nOrth',
        'CV.k',
        'CV.type',
        'tssx',
        'tssy'
      ),
      Value = c(center, scale, 1, nc - 1, cv.k, cv.type, tssx, tssy)
    )


    # define slots for OPLS_Torben object
    mod_opls = new('OPLS_MetaboMate',
      type = type,
      t_pred = pls_comp_all$scores,
      p_pred = pls_comp_all$loadings,
      w_pred = pls_comp_all$weights,
      betas_pred = pls_comp_all$betas,
      Qpc = pls_comp_all$Qpc,
      t_cv = cbind(t_cv[, (nc - 1)]),
      t_orth_cv = cbind(t_orth_cv[, (nc - 1)]) ,
      t_orth = cbind(output[['t_orth']][, 1:(nc - 1)]),
      p_orth = rbind(output[['p_orth']][1:(nc - 1), ]),
      w_orth = rbind(output[['w_orth']][1:(nc - 1), ]),
      nPC = nc - 1,
      summary = model.summary,
      X_orth = Xorth,
      Y_res = pls_comp_all$Y.res,
      Xcenter = attr(XcsTot, "scaled:center"),
      Xscale = attr(XcsTot, "scaled:scale"),
      Ycenter = attr(YcsTot, "scaled:center"),
      Yscale = attr(YcsTot, "scaled:scale"),
      Yout = Y_out[[2]],
      Parameters = dd,
      #sdFeatures=sdX,
      E = E
    )

    return(mod_opls)
  } else{
    print('No sign. orthogonal components - try PLS!')
    return(NULL)
  }

}
