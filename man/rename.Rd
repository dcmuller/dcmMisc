\name{rename}
\alias{rename}
\title{Rename columns of a data frame}
\usage{
rename(df, oldnames, newnames)
}
\arguments{
  \item{df}{a \code{data.frame} containing variables to be
  renamed}

  \item{oldnames}{character vector of one or more column
  names of \code{df}}

  \item{newnames}{character vector of new names for the
  columns listed in \code{oldnames}}
}
\value{
the \code{data.frame} df, with columns renamed
}
\description{
\code{rename} renames one or more columns of a data frame.
}

