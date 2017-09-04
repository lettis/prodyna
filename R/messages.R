#' Warning/Error Message.
#'
#' Get the warning/error message for the given type of warning/error.
#'
#' Possible warning/error types are:
#' \itemize{
#'   \item \code{caching} - Cached files exist. \cr
#'                          \code{arg} should give the function that produced
#'                          the warning.
#'   \item \code{noExec}  - Command could not be executed.
#' }
#'
#' @param type Character, type of warning/error
#' @param arg Additional information.
msg <- function(type, arg=NULL) {
  switch(type,
    caching = paste(arg, ": ",
                    "Output file(s) already exist. Nothing to compute.\n", sep=""),
    noExec  = paste("An error occured executing", arg, sep="\n"),
    format  = paste("The file", arg, "does not have the expected format."),
    missingFile = paste("File(s)", paste(arg, sep=", ", collapse = ", "), "do not exist."),
    "Unknown warning or error.\n"
    )
}
