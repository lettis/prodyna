
#' Radii estimation for density-based clustering
#' @param rc Reaction coordinates.
#' @param radii Vector of radii to test.
#' @export
clustering.estimate.radii <- function(rc, radii) {
  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
  if ( ! file.exists(rc)) {
    cls_dir <- get.fullPath(cls_dir)
  }
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
               "-p pop")
  message("running radii estimation for density-based clustering ...")
  system(cmd)
  message("... finished")
}



#' compute neighborhood for efficient density-based geometric clustering
#' @param rc Reaction coordinates to cluster.
#' @param radius Radius for density estimation.
#' @export
clustering.compute.neighborhood <- function(rc, radius) {
  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
  if ( ! file.exists(rc)) {
    cls_dir <- get.fullPath(cls_dir)
  }
  # run neighborhood computation
  cmd <- paste("cd", cls_dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-r",
               radius,
               paste("-D fe_", format(radius, nsmall=6), sep=""),
               "-b nn")
  message("running neighborhood computation for density-based clustering ...")
  system(cmd)
  message("... finished")
}



#' run screening (density-based geometric clustering)
#' @param rc Reaction coordinates to cluster.
#' @param radius Radius for density estimation.
#' @export
clustering.screening <- function(rc, radius) { #, min=NULL, max=NULL, step=NULL) {

  ##TODO handle min, max, step

  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
  if ( ! file.exists(rc)) {
    cls_dir <- get.fullPath(cls_dir)
  }
  # run neighborhood computation
  cmd <- paste("cd", cls_dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-r",
               radius,
               paste("-D fe_", format(radius, nsmall=6), sep=""),
               "-B nn",
               "-T -1",
               "-o clust")
  message("running screening ...")
  system(cmd)
  message("... finished")
}

#' construct density network from landscape screening
#' @param rc Reaction coordinates to cluster.
#' @param minpop Minimum population to count microstate as node.
#' @export
clustering.densityNetwork <- function(rc, minpop) {

  ## TODO: min, max, step

  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
  if ( ! file.exists(rc)) {
    cls_dir <- get.fullPath(cls_dir)
  }
  # run neighborhood computation
  cmd <- paste("cd", cls_dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "network",
               "-p",
               minpop)
  message("constructing density network ...")
  system(cmd)
  message("... finished")
}



#' construct density network from landscape screening
#' @param rc Reaction coordinates to cluster.
#' @param radius Density estimation radius.
#' @param sorted Order microstate names by descending populations. (default: true)
#' @param output Filename of resulting microstate trajectory. default: "microstates"
#' @export
clustering.microstates <- function(rc, radius, sorted=TRUE, output="microstates") {
  .check.projectPath()
  # create clustering directory
  suppressWarnings(cls_dir <- normalizePath(paste(rc, "clustering", sep=".")))
  if ( ! file.exists(rc)) {
    cls_dir <- get.fullPath(cls_dir)
  }
  # run neighborhood computation
  cmd <- paste("cd", cls_dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-r",
               radius,
               paste("-D fe_", format(radius, nsmall=6), sep=""),
               "-B nn",
               "-i network_end_node_traj.dat",
               paste("-o", output))
  message("computing microstates from density network ...")
  system(cmd)
  if (sorted) {
    ## reorder state names by descending population
    require(dplyr, quietly=TRUE, warn.conflicts=FALSE)
    output <- paste(cls_dir, output, sep="/")
    traj <- data.table::fread(output,
                              verbose=FALSE,
                              showProgress=FALSE)
    states <- unique(traj$V1)
    pops <- sapply(states, function(s){
      sum(traj$V1 == s)
    })
    # match old to new names
    matching <- data.frame(states, pops) %>%
      arrange(desc(pops)) %>%
      mutate(id=1:length(states))
    # convert to simple dictionary for fast retrieval
    dict <- list()
    for (s in states) {
      dict[s] <- matching %>%
        filter(states==s) %>%
        select(id)
    }
    # convert trajectory
    traj <- sapply(traj$V1, function(s){
      dict[[s]]
    })
    data.table::fwrite(list(traj),
                       output,
                       col.names=FALSE,
                       verbose=FALSE,
                       showProgress=FALSE)
  }
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

  #TODO find plotting method with higher performance
  #     (not geom_point, perhaps 1d hist with bins=1% of n data points)

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

