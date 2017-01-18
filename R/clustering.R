
#' Radii estimation for density-based clustering
#' @param rc Reaction coordinates.
#' @param radii Vector of radii to test.
#' @export
clustering.estimate.radii <- function(rc, radii) {
  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
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

#' Return data frame with population per frame
#' @param rc The clustered reaction coordinates
#' @param radii Selection of radii. If NULL (default), get all available
#' @export
clustering.get.pops <- function(rc, radii=NULL) {
  if (is.null(radii)) {
    popfiles <- list.files(get.fullPath(paste(rc, "clustering", sep=".")),
                           pattern="pop_*",
                           full.names=TRUE)
  } else {
    popfiles <- sapply(radii, function(r) {
      get.fullPath(c(paste(rc, "clustering", sep="."),
                     paste("pop_", sprintf("%0.6f", r), sep="")))
    })
  }
  do.call(data.frame,
          lapply(popfiles, function(fname){
            pops <- data.frame(read.table(fname)[[1]])
            colnames(pops) <- tail(strsplit(fname, split="/")[[1]], n=1)
            pops
          }))
}


#' Plot per-frame populations for given radii.
#' @param rc Reaction coordinates used for clustering.
#' @param radii Radii selection. If NULL (default), plot all available.
#' @param logy Plot with logarithmic y-scale.
#' @export
clustering.plot.pops <- function(rc, radii=NULL, logy=TRUE) {
  suppressMessages(library(ggplot2))
  suppressMessages(library(reshape2))
  pops <- clustering.get.pops(rc, radii)
  for(i in colnames(pops)) {
    pops[[i]] <- sort(pops[[i]], decreasing=TRUE)
  }
  pops$id <- 1:nrow(pops)
  p <- ggplot(melt(pops, id.vars="id")) +
        geom_line(aes(x=id, y=value, color=variable)) +
        xlab("frames") +
        ylab("sorted populations")
  if (logy) {
    p <- p + scale_y_log10()
  }
  p
}
