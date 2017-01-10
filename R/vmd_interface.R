
#' run VMD instance
#'
#' Run VMD in background, starting a remote control server on localhost:5555.
#' @export
run.vmd <- function() {
  invisible(subprocess::spawn_process(get.binary("vmd"),
                                      arguments = c("-e",
                                                    system.file("extdata",
                                                                "remote_ctl.tcl",
                                                                package="prodyna"))))
}

#' exec Tcl command
#' @export
exec.tcl <- function(tcl) {
  con <- socketConnection(port=5555)
  write(tcl, con)
}


# TODO: more commands for direct manipulation.
#    use 'logfile console' to get Tcl commands for gui actions


