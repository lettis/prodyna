
.default_binaries <- list("gmx"="/usr/local/gromacs/bin/gmx",
                          "fastpca"="/usr/local/bin/fastpca",
                          "clustering"="/usr/local/bin/clustering",
                          "bash"="/bin/bash",
                          "awk"="/usr/bin/awk",
                          "paste"="/usr/bin/paste")

.binaries <- new.env()

#' Set path to binary
#'
#' @param name Name of the binary
#' @param path Absolute system path to binary
#' @examples
#' set.binary("bash", "/bin/bash")
#' @export
set.binary <- function(name, path) {
  assign(name, path, envir=.binaries)
}

#' Get path to binary
#'
#' @param name Name of the binary
#' @export
get.binary <- function(name) {
  if (exists(name, envir=.binaries, inherits=FALSE)) {
    get(name, envir=.binaries)
  } else {
    .default_binaries[[name]]
  }
}

#' List names and paths to binaries called by prodyna
#'
#' @export
list.binaries <- function() {
  lapply(names(.default_binaries), function(n) {
    c(n, get.binary(n))
  })
}

#' Project list
#'
#' Holds all project-related information.
#' Automatically updated by generator functions.
#'
#' @export
project <- list()



