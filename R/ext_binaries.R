
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
#' Execution is stopped if command gives a non-zero exit status.
#'
#' @param cmdName Character, name of binary to run.
#' @param args Character vector, arguments to the command.
run.cmd <- function(cmd_name, args) {
  tryCatch(
    system2(get.binary(cmd_name), args, stdout = TRUE, stderr = TRUE),
    error   = function(c) {stop(msg.err("noExec", cmd_name), call. = FALSE)},
    warning = function(c) {stop(c[[1]], call. = FALSE)}
  )
}

#' Run commands.
#'
#' Runs the string of commands given by cmds.
#' Execution is stopped if commands give a non-zero exit status.
#'
#' Use \code{run.cmd} for a single command.
#'
#' @param cmds  Character, name of binary to run.
run.cmds <- function(cmds) {
  tryCatch(
    system(cmds, intern = TRUE),
    error   = function(c) {stop(msg.err("noExec", cmds), call. = FALSE)},
    warning = function(c) {stop(c[[1]], call. = FALSE)}
    )
}
