
#' Warning messages.
#'
.warnings <- list(
  "no_init" = "No project has been initialized so far. Call prodyna::initialize(...)",
  "diverging_path" = "Current directory and project paths diverge. Resetting current directory to project path.",
  "missing_files" = "Necessary files are missing.")

#' Error messages.
#'
#' @param type Character, type of error.
#'             `noExec` Command could not be exectued.
#' @param arg  Additional information.
#' @return Error message.
msg.err <- function(type, arg=NULL) {
  switch(type,
         noExec = paste("Could not execute ", arg,
                        ". Adapt path to binary using 'set.binary'.", sep=""),
         "Unknown error occured.")
}
