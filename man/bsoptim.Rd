\name{bsoptim}
\alias{bsoptim}
\title{Bootstrap optimism corrected c-statistic for lasso variable selection}
\usage{
bsoptim(fit, reps)
}
\arguments{
  \item{fit}{object returned by \code{lassofit}}

  \item{reps}{number of bootstrap replications}
}
\value{
a numeric vector containing the apparent c-statistic, its
estimated optimism, and the corrected c-statistic
}
\description{
\code{bsoptim} takes an object produced by \code{lassofit}
and calculates bootstrap optimism corrected model
discrimination.
}
\author{
David C Muller
}

