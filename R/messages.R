
#' Warning messages.
#'
.warnings <- list(
  "no_init" = "No project has been initialized so far. Call prodyna::initialize(...)",
  "diverging_path" = "Current directory and project paths diverge. Resetting current directory to project path.",
  "missing_files" = "Necessary files are missing.")

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
    noExec  = paste("Could not execute ", arg,
                    ". Adapt path to binary using 'set.binary'.\n", sep=""),
    "Unknown warning or error.\n"
    )
}
