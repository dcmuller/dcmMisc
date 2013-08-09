#' @title Rename columns of a data frame
#' 
#' @description \code{rename} renames one or more columns of a data frame.
#' 
#' @param df a \code{data.frame} containing variables to be renamed
#' @param oldnames character vector of one or more column names of \code{df}
#' @param newnames character vector of new names for the columns listed in
#'   \code{oldnames}
#'   
#' @return the \code{data.frame} df, with columns renamed
#' @export
rename <- function(df, oldnames, newnames) {
  if (length(oldnames) != length(newnames)) {
    stop("oldnames must be the same length as newnames")
  }
  if (!is.data.frame(df)) {
    stop("First argument is not a data.frame")
  }
for (i in 1:length(oldnames))
  names(df)[names(df) == oldnames[i]] <- newnames[i]
return(df)
}