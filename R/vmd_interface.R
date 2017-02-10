

# TODO: more commands for direct manipulation.
#    use 'logfile console' to get Tcl commands for gui actions


#' run VMD instance
#'
#' Run VMD in background, starting a remote control server on localhost:5555.
#' @export
vmd.run <- function() {
  invisible(subprocess::spawn_process(get.binary("vmd"),
                                      arguments = c("-e",
                                                    system.file("extdata",
                                                                "remote_ctl.tcl",
                                                                package="prodyna"))))
}

#' exec VMD-Tcl command
#' @param tcl The TcL command as string.
#' @export
vmd.exec <- function(tcl) {
  con <- socketConnection(port=5555)
  write(tcl, con)
  write("close", con)
  Sys.sleep(1)
  close(con)
}

#' close the running VMD instance
#' @export
vmd.quit <- function() {
  vmd.exec("quit")
}

#' load pdb file into VMD
#' @param filename The pdb file.
#' @export
vmd.loadPdb <- function(filename) {
  vmd.exec(paste("mol new {", filename, "}", sep=''))
}

#' set rendering style
#'
#' @param style The rendering style. One of: "Licorice"
#' @export
vmd.setStyle <- function(style, molid=0) {
  style_opts <- list()
  style_opts["Licorice"] <-"0.300000 12.000000 12.000000"
  vmd.exec(paste("mol modstyle 0",
                 molid,
                 style,
                 style_opts[[style]]))
}

#' render current view to postscript
#' @param filename Filename of rendered output
#' @export
vmd.renderPostscript <- function(filename="vmdscene.ps") {
  vmd.exec("color Display Background white")
  vmd.exec(paste("render PostScript", filename))
}
