# ------------------------------------------------------------------------
# file: 00_project_setup.R author: Jason Grafmiller date:
# 2017/07/24

# Description: R code for conducting random forest analysis
# as reported in Grafmiller, J & B. Szmrecsanyi.
# 'Mapping out particle placement in Englishes around the
# world'. Language Variation and Change.

# ------------------------------------------------------------------------


# load core libraries -----------------------------------------------------

# for data manipulation
library(tidyverse)
library(reshape2)

# for plotting
library(gridExtra)
library(scales)
library(RColorBrewer)

# for regression modelling, the following packages should be
# installed: - lme4, rms, Hmisc, MuMIn, merTools, car

# for random forest modelling, the following packages should
# be installed: - party, edarf

# for neighborNet modelling, the following packages should be
# installed: - phangorn


# custom functions --------------------------------------------------------

# Custom functions

notify <- function(x = NULL) {
  # Generate a pop-up notification window for use with processes that
  # take long time to run, e.g. varimp() with random forests
  if (is.null(x)) {
    system("CMD /C \"ECHO The R process has finished running && PAUSE\"",
      invisible = FALSE, wait = FALSE)
  } else {
    system(paste("CMD /C \"ECHO ", x, " && PAUSE\""), invisible = FALSE,
      wait = FALSE)
  }
}

# centering and standardizing functions following Gelman (2008)
c. <- function(x) {
  y <- as.numeric(x)
  return (y - mean(y))
}

z. <- function(x, factor = 1) {
  v <- (x - mean(x))/(factor*sd(x))
  return (v)
}

# --- For working with merMod objects:

collin.fnc.mer <- function(fit, ...) {
  # adaption of languageR::collin.fnc() for compatibility with
  # current version of lme4 merMod objects and R 3.4
  require(languageR, quietly = T)
  if (class(fit) == "lmerMod" || class(fit) == "glmerMod") {
    data <- getME(fit, "X")[,-1]
    colvector <- 1:ncol(data)
    std.fnc <- function(vec) (vec - mean(vec))/sqrt(var(vec))
    # New from R 3.4: Add as.vector() to avoid warning 'Recycling
    # array of length 1 in vector array arithmetic is deprecated'
    scale.fnc <- function(vec) (vec/sqrt(as.vector(t(vec) %*%
      vec)))
    x = data[, colvector]
    xlengte = length(x[, 1])
    colnames = dimnames(x)[[2]]
    onevec = rep(1, xlengte)
    Xdesign = cbind(onevec, as.matrix(x))
    X = Xdesign
    ncols = length(X[1, ])
    for (i in 1:ncols) {
      X[, i] = scale.fnc(as.numeric(X[, i]))
    }
    svd.X = svd(X, nu = 0)
    nu.X = max(svd.X$d)/svd.X$d
    kappa.X = max(nu.X)
    pi.X = svd.X$v
    for (i in (1:length(svd.X$d))) {
      pi.X[, i] = svd.X$v[, i]/svd.X$d[i]
      pi.X[, i] = pi.X[, i]^2
    }
    for (i in 1:length(svd.X$d)) {
      pi.X[i, ] = pi.X[i, ]/sum(pi.X[i, ])
    }
    pi.X = t(pi.X)
    pi.X = as.data.frame(pi.X)
    dimnames(pi.X)[[2]] = c("Constant", colnames)
    return (list(svd = svd.X, cindex = nu.X, cnumber = kappa.X,
      pi = pi.X))
  } else stop("model not a merMod object")
}

overdisp_fun <- function(fit, ...) {
  # Diagnose overdispersion in a model's response variable
  # number of variance parameters in an n-by-n
  # variance-covariance matrix. Adapted from:
  # http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  model.df <- sum(sapply(VarCorr(fit), vpars)) + length(fixef(fit))
  rdf <- nrow(model.frame(fit)) - model.df
  rp <- residuals(fit, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  data.frame(
    chisq = Pearson.chisq,
    ratio = prat,
    rdf = rdf,
    p = pval)
}

scores_mer <- function(fit, R2 = F, digits = 3, ...) {
  # Computes a number of model accuracy statistics:
  # - N = number of observations
  # - C = concordance index
  # - Dxy = Somer's Dxy
  # - AICc = corrected AIC (see Burnham & Anderson 2002)
  # - kappa = measure of data multicollinearity
  # - predicted.corr = Percentage of observations predicted correctly by
  #   the model
  # - baseline = Baseline accuracy (proportion of most common response)
  # - Log.Loss = (average log-likelihood)
  # - Class.1.acc = Percentage of observations of response 1 that were
  #   predicted correctly by the model (here the continuous PV variant)
  # - Class.2.acc = Percentage of observations of response 2 that were
  #   predicted correctly by the model (here the split PV variant)
  # - Avg.per.Class = Mean per class accuracy (mean of Class.1.acc and
  #   Class.2.acc)
  # - R2.m = Marginal R2
  # - R2.c = Conditional R2
  require(lme4)
  require(MuMIn)
  d <- fit@frame
  probs <- predict(fit, type = "response")
  preds <- ifelse(probs > 0.5, 1, 0)
  y <- getME(fit, "y")
  d2 <- cbind(d, data.frame(probs = probs, preds = preds, y = y))
  # AUC and Dxy
  C <- Hmisc::somers2(probs, y)[1] %>% as.vector
  Dxy <- Hmisc::somers2(probs, y)[2] %>% as.vector
  AICc <- MuMIn::AICc(fit)
  logL <- logLik(fit)[1]
  N <- nrow(d)
  kappa <- collin.fnc.mer(fit)$cnumber
  # log-loss (- average log-likelihood)
  LL <- -mean(y * log(probs) + (1 - y) * log(1 - probs))
  # predictive accuracies
  acc <- length(which(preds == y))/nrow(d)
  base <- 1 - mean(y)
  resp <- levels(d[, 1])
  resp_col <- names(d)[1]
  d.1 <- dplyr::filter_(d2, paste0(resp_col, " == '", resp[1],
    "'")) %>% droplevels
  d.2 <- dplyr::filter_(d2, paste0(resp_col, " == '", resp[2],
    "'")) %>% droplevels
  class1.acc <- length(which(d.1$preds == d.1$y))/nrow(d.1)
  class2.acc <- length(which(d.2$preds == d.2$y))/nrow(d.2)
  out <- data.frame(N = N, C = C, Dxy = Dxy, AICc = AICc, kappa = kappa,
    predicted.corr = acc, baseline = base, Log.Loss = LL,
    Class.1.acc = class1.acc, Class.2.acc = class2.acc,
    Avg.Per.Class = mean(c(class1.acc, class2.acc)))
  if (R2) {
    R2 <- MuMIn::r.squaredGLMM(fit)
    out <- mutate(out, R2.m = R2[[1]], R2.c = R2[[2]])
  }
  if (is.null(digits)) {
    return(out)
  } else return (round(out, digits))
}


VIF_mer <- function(fit) {
  # Compute variance inflation factor for individual predictors
  # in a lme4 model. Adapted from rms::vif()
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v2 <- diag(solve(v/(d %o% d)))
  names(v2) <- nam
  v2 <- sort(v2, decreasing = TRUE)
  return(v2)
}

permute_varimp_mer <- function(fit, verbose = FALSE, ranef = FALSE){
  # Calculates the differences in model fit scores between a full model and
  # models with each individual predictor permuted. For each predictor in the
  # model, the values of that predictor are randomly permuted to break their
  # association with the response, and the model is re-fit to a new dataset
  # containing the premuted values. The fit of the new model is compared to
  # that of the original model. See Baayen (2011:308).
  # Three measures of model fit are calculated: The concordance C,
  # Accuracy (% correctly predicted), and AICc.
  require(lme4, quietly = T)
  require(Hmisc, quietly = T)
  require(MuMIn, quietly = T)
  # Want this to work with glm and glmer models...
  if (class(fit)[1] %in% c("lmerMod", "glmerMod")){
    y <- fit@resp$y
    data <- as.data.frame(fit@frame)
    if (ranef){
      random <- names(ranef(fit))
      vars <- colnames(attr(attr(fit@frame, "terms"), "factors"))
    } else {
      # only fixed effs
      vars <- strsplit(toString(attr(attr(fit@frame, "terms"),
       "predvars.fixed")), ', ')[[1]][-c(1:2)]
    }
  } else if (class(fit)[1] == "lrm"){
    vars <- colnames(attr(f.lrm$terms, "factors"))
    if(is.null(fit$y)) {
      stop("Must specify lrm(..., y = T) when fitting lrm().")
      } else {
        y <- as.numeric(fit$y) - 1
      }
    ranef <- FALSE
  } else if (class(fit)[1] == "glm"){
    data <- fit$data
    vars <- colnames(attr(attr(fit$model, "terms"), "factors"))
    y <- as.numeric(fit$y)
    ranef <- FALSE
  } else {
    stop(paste(fit, "is not a supported class"))
  }
  # remove interaction terms
  vars <- vars[grep(":", vars, invert = T)]
  # get predictions from full model
  full.probs <- predict(fit, type = "response")
  full.C <- somers2(full.probs, y)[[1]]
  full.acc <- mean(round(full.probs) == y)
  full.AICc <- MuMIn::AICc(fit)
  varimp_mat <- matrix(nrow = length(vars), ncol = 3)
  if (verbose) {cat("variables run:\n")}
  # loop through (fixed effects) predictors
  for (i in seq(1, length(vars))){
    # find main effect and any interactions
    d <- data
    d[, vars[i]] <- sample(data[, vars[i]]) # reshuffle values
    if (class(fit)[1] == "glmerMod"){
      # Make sure the updated model inherits the control settings from original
      ctrl = glmerControl(optimizer = fit@optinfo$optimizer,
                          optCtrl = fit@optinfo$control)
      new_fit <- update(fit, data = d, control = ctrl)
    } else if (class(fit)[1] == "lmerMod") {
      # Make sure the updated model inherits the control settings from original
      ctrl = lmerControl(optimizer = fit@optinfo$optimizer,
                         optCtrl = fit@optinfo$control)
      new_fit <- update(fit, data = d, control = ctrl)
    } else {
      new_fit <- update(fit, data = d)
    }
    new.probs <- predict(new_fit, type = "response")

    if (class(fit)[1] %in% c("lmerMod", "glmerMod")){
      new_y <- new_fit@resp$y
    } else {
      new_y <- new_fit$y
    }
    new.C <- somers2(new.probs, new_y)[[1]]
    C.diff <- full.C - new.C
    new.acc <- mean(round(new.probs) == new_y)
    Acc.diff <- full.acc - new.acc
    new.AICc <- MuMIn::AICc(new_fit)
    # The difference in AICc here is the same as the likelihood ratio
    AICc.diff <- new.AICc - full.AICc
    varimp_mat[i, ] <- c(C.diff, Acc.diff, AICc.diff)
    if (verbose) {cat(vars[i], "... ", sep = "")}
  }
  rownames(varimp_mat) <- vars
  colnames(varimp_mat) <- c("C", "accuracy", "AICc")
  return(as.data.frame(varimp_mat))
}


# --- For working with ggplot

grid.share.legend <- function(..., ncol = length(list(...)),
  nrow = 1, position = c("bottom", "right")) {
  # Allow multiple plots to be printed in a single window with a shared legend.
  # See:
  # https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  require(grid)
  require(gridExtra)
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position, bottom = arrangeGrob(do.call(arrangeGrob,
    gl), legend, ncol = 1, heights = unit.c(unit(1, "npc") -
    lheight, lheight)), right = arrangeGrob(do.call(arrangeGrob,
    gl), legend, ncol = 2, widths = unit.c(unit(1, "npc") -
    lwidth, lwidth)))
  # grid.newpage()
  grid.draw(combined)
}


# references --------------------------------------------------------------

# Baayen, R. Harald. 2011. Corpus linguistics and Naive Discriminative Learning.
#   Revista Brasileira de Linguística Aplicada 11(2). 295–328.

# Burnham, Kenneth P. & David R. Anderson. 2002. Model Selection and Multimodel
#   Inference: A Practical Information-Theoretic Approach. New York: Springer.

# Gelman, Andrew. 2008. Scaling Regression Inputs by Dividing by Two Standard
#   Deviations. Statistics in Medicine 27(15). 2865–2873. doi:10.1002/sim.3107.

