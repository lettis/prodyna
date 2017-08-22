
.default_binaries <- list("gmx"="/usr/local/gromacs/bin/gmx",
                          "vmd"="/usr/local/bin/vmd",
                          "fastpca"="/usr/local/bin/fastpca",
                          "clustering"="/usr/local/bin/clustering",
                          "entropical"="/usr/local/bin/entropical",
                          "bash"="/bin/bash",
                          "awk"="/usr/bin/awk",
                          "paste"="/usr/bin/paste")

.binaries <- new.env()

#' Set path to binary.
#'
#' @param name Name of the binary.
#' @param path Absolute system path to binary.
#' @examples
#' set.binary("bash", "/bin/bash")
#' @export
set.binary <- function(name, path) {
  assign(name, path, envir=.binaries)
}

#' Get path to binary.
#'
#' @param name Name of the binary.
#' @export
get.binary <- function(name) {
  if (exists(name, envir=.binaries, inherits=FALSE)) {
    get(name, envir=.binaries)
  } else {
    .default_binaries[[name]]
  }
}

#' List of external binaries.
#'
#' List names and paths to binaries called by prodyna.
#' @export
list.binaries <- function() {
  lapply(names(.default_binaries), function(n) {
    c(n, get.binary(n))
  })
}

#' Run single command.
#'
#' Runs command with given arguments.
#'
#' Wrapper for \code{\link{system2}} that stops the execution and throws an
#' error if the command gives a non-zero exit status.
#'
#' @param cmd Character, command to run.
#' @param args Character vector, arguments to the command.
#' @param onError Function, function to be performed in case of an error.
run.cmd <- function(cmd, args, onError=NULL) {
  errorAction <- function(c) {
    if (!is.null(onError)) {
      onError()
    }
    stop(msg("noExec", paste(cmd, paste(args, collapse =" "))),
         call. = F)
  }
  tryCatch(
    system2(cmd, args, stdout = TRUE, stderr = TRUE),
    error   = errorAction,
    warning = errorAction
  )
}

#' Run commands.
#'
#' Runs the string of commands given by cmds.
#'
#' Wrapper for \code{\link{system}} that stops the execution and throws an error
#' if the command gives a non-zero exit status.
#' Use \code{\link{run.cmd}} for a single command.
#'
#' @param cmds  Character, string of commands to be run.
#' @param onError Function, function to be performed in case of an error.
run.cmds <- function(cmds, onError=NULL) {
  errorAction <- function(c) {
    if (!is.null(onError)) {
      onError()
    }
    stop(msg("noExec", cmds), call. = F)
  }
  output = tryCatch(
    expr    = system(cmds, intern = TRUE),
    error   = errorAction,
    warning = errorAction
    )
}
