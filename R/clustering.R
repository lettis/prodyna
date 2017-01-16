
#' Radii estimation for density-based clustering
#' @param rc Reaction coordinates.
#' @param radii Vector of radii to test.
#' @export
clustering.estimate.radii <- function(rc, radii) {
  .check.projectPath()
  # create clustering directory
  cls_dir <- normalizePath(paste(rc, "clustering", sep="."))
  if ( ! dir.exists(cls_dir)) {
    dir.create(cls_dir)
  }
  # create link to reaction coords
  rc_link <- paste(cls_dir, "reaction_coords", sep="/")
  if ( ! file.exists(rc_link)) {
    file.symlink(paste("../", rc, sep=""), rc_link)
  }
  # run clustering
  cmd <- paste("cd", cls_dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-R",
               paste(radii, collapse = " "),
               "-d fe",
               "-p pop",
               "-v")
  message("running radii estimation for density-based clustering ...")
  system(cmd)
  message("... finished")
}


#' Plot per-frame populations for given radii.
#' @param rc Reaction coordinates used for clustering.
#' @param radii Radii selection. If NULL (default), plot all available.
#' @param logy Plot with logarithmic y-scale.
clustering.plot.pops <- function(rc, radii=NULL, logy=TRUE) {
  #TODO implement
}
