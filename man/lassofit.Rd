\name{lassofit}
\alias{lassofit}
\title{Variable selection via lasso}
\usage{
lassofit(x, y, family, nfolds = 10)
}
\arguments{
  \item{x}{full design matrix}

  \item{y}{vector of responses ("dependent variable")}

  \item{family}{outcome family (see \code{\link{glm}})}

  \item{nfolds}{number of folds for cross validation of
  penalty term (see \code{\link{glmnet}})}
}
\value{
unconstrained glm of y on the selected subset of columns of
x
}
\description{
\code{lassofit} uses \code{glmnet} to do variable selection
via the lasso, and then returns an unconstrained glm which
includes the selected variables. The returned object can be
passed to \code{bsoptim} which calculates bootstrap
optimism corrected indices of discrimination.
}
\author{
David C Muller
}

