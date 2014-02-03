#' @title Variable selection via lasso
#' 
#' @description \code{lassofit} uses \code{glmnet} to do variable selection via
#' the lasso, and then returns an unconstrained glm which includes the selected
#' variables. The returned object can be passed to \code{bsoptim} which 
#' calculates bootstrap optimism corrected indices of discrimination.
#' 
#' @param x full design matrix
#' @param y vector of responses ("dependent variable") 
#' @param family outcome family (see \code{\link{glm}})
#' @param nfolds number of folds for cross validation of penalty term 
#'    (see \code{\link{glmnet}})
#'   
#' @return unconstrained glm of y on the selected subset of columns of x
#'
#' @author David C Muller
#' 
#' @import glmnet
#' @import rms
#' @export
lassofit <- function(x, y, family, nfolds=10) {
  call <- match.call()
  ## get column names from x
  if (is.null(colnames(x))) {
    nms <- paste0("X", 1:ncol(x))
    colnames(x) <- nms
  }
  xnames <- colnames(x)
  
  ## cross validated lasso to select shrinkage parameter (lambda)
  cv_fit <- cv.glmnet(x=x, y=y, family=family, nfolds=nfolds)
  ## get reduced fit corresponding to lambda.min
  full_fit <- glmnet(x=x, y=y, family=family)
  reduced_fit <- coef(full_fit, s=cv_fit$lambda.min)
  active_indx <- which(as.logical(reduced_fit != 0))
  active_vars <- rownames(reduced_fit)[active_indx][-1]
  
  ## if there are no active variables, fit an intercept only model
  if (length(active_vars)==0)
    active_vars <- "1"
  
  ## set up data.frame and formula for "unbiased" logistic fit
  ## using the selected variables
  dframe <- data.frame(y, x)
  names(dframe) <- c("y", xnames)
  m_formula <- formula(paste("y ~ ", paste(active_vars, collapse=" + ")))
  m_fit <- glm(m_formula, family=family, x=TRUE, y=TRUE, model=TRUE, data = dframe)
  m_fit$call <- call
  m_fit$lassofit_x <- x
  m_fit$lassofit_y <- y
  m_fit$lassofit_model <- dframe
  m_fit  
}

#' @title Bootstrap optimism corrected c-statistic for lasso variable selection
#' 
#' @description \code{bsoptim} takes an object produced by \code{lassofit} and
#' calculates bootstrap optimism corrected model discrimination. 
#' 
#' @param fit object returned by \code{lassofit}
#' @param reps number of bootstrap replications
#'   
#' @return a numeric vector containing the apparent c-statistic, its estimated
#' optimism, and the corrected c-statistic
#'
#' @author David C Muller
#' 
#' @import glmnet
#' @import rms
#' @export
bsoptim <- function(fit, reps) {
  ## initialise optimism to 0
  o <- 0
  sampsize <- nrow(fit$model)
  call <- fit$call
  
  ## apparent c from full model
  c_app <- somers2(x=predict(fit, type="response"), fit$lassofit_y)[1]
  
  ## loop over bootstrap samples
  i <- 1
  while (i <= reps) {
    bsindex <- sample(1:sampsize, sampsize, replace=TRUE)
    ynew <- fit$lassofit_y[bsindex]
    xnew <- fit$lassofit_x[bsindex, ]
    bsfit <- update(fit, x=xnew, y=ynew)
    ## apparent c from the bootstrap sample
    c_boot <- somers2(x=predict(bsfit, type="response"), bsfit$y)[1]
    
    ## out of sample prediction
    c_orig <- somers2(x=predict(bsfit, newdata=fit$lassofit_model, type="response"),
                      y=fit$y)[1]
    o <- o + (c_boot - c_orig)
    i <- i + 1
  }
  o <- o/reps
  corrected <- c_app - o
  res <- c(corrected, o, c_app, reps)
  names(res) <- c("corrected", "optimism", "apparent", "N_reps")
  res
}


## Example (not run)
# set.seed(137494)
# y <-rnorm(200)
# x1 <- matrix(rnorm(200*20),200,20)
# x2 <- matrix(y+5*rnorm(200*10),200,10)
# x <- cbind(x1,x2)
# y <- as.numeric(y>0)
# 
# myfit <- lassofit(x=x, y=y, family="binomial")
# my_validated <- bsoptim(fit=myfit, reps=50)
# print(my_validated)
