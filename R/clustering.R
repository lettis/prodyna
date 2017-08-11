
#' Radii estimation for density-based clustering.
#'
#' @param coords Character, name of the reaction coordinates file.
#' @param radii Numeric vector, radii to test.
#' @param dir Character, name of the output directory. If \code{NULL} the output
#'  directory is <coords>.clustering. \cr
#'  If the ouput directory is created if it does not exist already.
#' @export
clustering.estimate.radii <- function(coords, radii, dir=NULL) {

  if (is.null(dir)) {
    dir <- normalizePath(paste(coords, "clustering", sep="."))
  }
  # create clustering directory
  if (!dir.exists(dir)) {
    message(paste("Creating output directory", dir))
    dir.create(dir)
  }
  # create link to reaction coords
  coords_link <- paste(dir, "reaction_coords", sep="/")
  if (!file.exists(coords_link)) {
    file.symlink(normalizePath(coords), coords_link)
  }
  # run clustering
  cmd <- paste("cd", dir, ";")
  cmd <- paste(cmd,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-R",
               paste(radii, collapse = " "),
               "-d fe",
               "-p pop")
  message("running radii estimation for density-based clustering ...")
  run.cmds(cmd)
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

#' rename microstates from 1:N, sorted by descending populations
#' @param microstates Either filename or vector of microstate trajectory
#' @param output Filename of renamed microstate trajectory.
#'               (default: NULL, i.e. just return microstate trajectory as vector)
#' @export
clustering.microstates.renamed <- function(microstates, output=NULL) {
  require(dplyr)
  if (is.character(microstates)) {
    microstates <- data.table::fread(microstates,
                                     verbose=FALSE,
                                     showProgress=FALSE)[[1]]
  }
  microstates_renamed <- microstates
  n_microstates <- length(unique(microstates))
  counts <- data.frame(state=microstates) %>%
    count(state) %>%
    arrange(desc(n)) %>%
    mutate(idx=1:n_microstates)
  for (state in unique(microstates)) {
    microstates_renamed[microstates==state] <- counts[counts$state==state,]$idx
  }
  if ( ! is.null(output)) {
    data.table::fwrite(data.frame(microstates_renamed),
                       sep=" ",
                       col.names=FALSE,
                       verbose=FALSE,
                       showProgress=FALSE)
  }

  microstates_renamed
}

#' Return data frame with population per frame
#' @param rc Either the clustered reaction coordinates or a path to the pop files
#' @param radii Selection of radii. If NULL (default), get all available
#' @export
clustering.get.pops <- function(rc, radii=NULL) {
  if (dir.exists(rc)) {
    rc_dir <- rc
  } else {
    rc_dir <- get.fullPath(paste(rc, "clustering", sep="."))
  }
  if (is.null(radii)) {
    popfiles <- list.files(rc_dir,
                           pattern="pop_*",
                           full.names=TRUE)
  } else {
    popfiles <- sapply(radii, function(r) {
      get.fullPath(c(rc_dir,
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
#' @param rc Either the reaction coordinates used for clustering or
#'           a path to the pop files.
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


#' Plot hierarchical network of MPP lumping
#' @param dirname Directory of MPP run
#' @export
clustering.plot.mppNetwork <- function(dirname) {
  require(igraph, quietly=TRUE, warn.conflicts=FALSE)
  require(dplyr, quietly=TRUE, warn.conflicts=FALSE)

  #### get the data
  # helper to get different data files from MPP directory
  get_mpp_data <- function(fname) {
    data.table::fread(paste(dirname,
                            fname,
                            sep="/"),
                      verbose=FALSE,
                      showProgress=FALSE)
  }
  # get populations at different qmin levels
  mpp_pops <- lapply(list.files(dirname, pattern="mpp_pop*"), function(fname){
    pops <- get_mpp_data(fname)
    colnames(pops) <- c("state", "pop")
    pops
  })
  # helper function to get max pops from union of p1 and p2
  get_max_pops <- function(p1, p2) {
    all_states <- unique(c(p1$state, p2$state))
    pops <- c()
    for (s in all_states) {
      pops <- c(pops, max(p1$pop[p1$state==s],
                          p2$pop[p2$state==s],
                          na.rm=TRUE))
    }
    data.frame(state=all_states, pop=pops)
  }
  # get total max pops for all states (i.e. pops before lumping)
  max_pops <- Reduce(get_max_pops, mpp_pops)
  # get transitions (i.e. at what qmin level are states lumped?)
  transitions <- get_mpp_data("mpp_transitions.dat")
  colnames(transitions) <- c("from", "to", "qmin")
  transitions$qmin <- round(transitions$qmin, digits=4)

  #### construct the network
  pops <- (max_pops %>% arrange(state) %>% select(pop))[[1]]
  g <- graph(edges=as.vector(rbind(transitions$from, transitions$to)))
  g$weights <- 5*transitions$qmin
  plot(g,
       vertex.color="white",
       vertex.size=log(pops),
       edge.arrow.size=0.2,
       edge.width=g$weights,
       layout=layout_with_kk)
}

