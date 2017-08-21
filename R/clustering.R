#' Prepare the clustering output directory.
#'
#' The clustering output directory is created if it does not exist already.
#' A file named 'reaction_coords' that is linked to the reaction coordinates
#' file is created within that directoy.
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#' @return absolute path to the output directory

clustering.prepareOutputDir <- function(rc, dir) {
  rc <- normalizePath(rc, mustWork=TRUE)

  if (is.null(dir)) {
    dir <- paste(rc, "clustering", sep=".")
  }
  # create clustering directory
  if (!dir.exists(dir)) {
    message(paste("Creating output directory", dir))
    dir.create(dir)
  }
  # create link to reaction coords
  coords_link <- paste(dir, "reaction_coords", sep="/")

  if (file.exists(coords_link)) {
    # check whether coords_link is linked to rc
    if(rc != Sys.readlink(coords_link)) {
      stop(paste("File reaction_coords already exists in the output directory ",
                 dir,
                 ", but is not linked to the given reaction coordinates file",
                 rc,
                 sep=""))
    }
  } else {
    file.symlink(rc, coords_link)
  }
  return(dir)
}


#' Radii estimation for density-based clustering.
#'
#' TODO: description
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radii Numeric vector, radii to test.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @export
clustering.estimate.radii <- function(rc, radii, dir=NULL) {

  dir <- clustering.prepareOutputDir(rc, dir)

  # run clustering
  cmds <- paste("cd", dir, ";")
  cmds <- paste(cmds,
                get.binary("clustering"),
                "density",
                "-f reaction_coords",
                "-R",
                paste(radii, collapse = " "),
                "-d fe",
                "-p pop")

  message("Running radii estimation for density-based clustering ...")

  run.cmds(cmds)

  message("... finished.")
}


#' Neighbourhood computation for efficient density-based geometric clustering.
#'
#' TODO: description
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radius Numeric, radius for density estimation.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <coords>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @export
clustering.compute.neighborhood <- function(rc, radius, dir=NULL) {

  dir <- clustering.prepareOutputDir(rc, dir)

  # run neighborhood computation
  cmds <- paste("cd", dir, ";")
  cmds <- paste(cmds,
                get.binary("clustering"),
                "density",
                "-f reaction_coords",
                "-r",
                radius,
                paste("-D fe_", format(radius, nsmall=6), sep=""),
                "-b nn")

  message("Running neighborhood computation for density-based clustering ...")

  run.cmds(cmds)

  message("... finished.")
}


#' Screening (density-based geometric clustering).
#'
#' TODO: description
#'
#' @param rc Character, name of the teaction coordinates file.
#' @param radius Numeric, radius for density estimation.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @param min
#' @param max
#' @param step
#' @export
clustering.screening <- function(rc, radius, dir=NULL, min=NULL, max=NULL, step=NULL) {

  ##TODO handle min, max, step

  dir <- clustering.prepareOutputDir(rc, dir)

  # run neighborhood computation
  cmds <- paste("cd", dir, ";")
  cmds <- paste(cmds,
                get.binary("clustering"),
                "density",
                "-f reaction_coords",
                "-r",
                radius,
                paste("-D fe_", format(radius, nsmall=6), sep=""),
                "-B nn",
                "-T -1",
                "-o clust")

  message("Running screening ...")

  run.cmds(cmds)

  message("... finished.")
}


#' Construct density network from landscape screening.
#'
#' TODO: description
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param minpop Numeric, minimum population to count microstate as node.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @param min
#' @param max
#' @param step
#' @export
clustering.densityNetwork <- function(rc, minpop, dir=NULL, min=NULL, max=NULL, step=NULL) {

  ## TODO: min, max, step

  dir <- clustering.prepareOutputDir(rc, dir)

  # construct density network
  cmds <- paste("cd", dir, ";")
  cmds <- paste(cmds,
                get.binary("clustering"),
                "network",
                "-p",
                minpop)

  message("Constructing density network ...")

  run.cmds(cmds)

  message("... finished.")
}


#' Compute Microstate trajectory.
#'
#' TODO: description
#'
#' @param rc Character, name of reaction coordinates file.
#' @param radius Numeric, density estimation radius.
#' @param sorted Logical, should microstate names be ordered by descending
#'  populations? (Default: \code{TRUE})
#' @param output Character, name of the resulting microstate trajectory file.
#'  (Default: "microstates")
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @importFrom dplyr desc arrange mutate filter select
#' @importFrom magrittr "%>%"
#' @export
clustering.microstates <- function(rc, radius, sorted=TRUE, output="microstates", dir=NULL) {
  dir <- clustering.prepareOutputDir(rc, dir)

  # run neighborhood computation
  cmds <- paste("cd", dir, ";")
  cmds <- paste(cmds,
               get.binary("clustering"),
               "density",
               "-f reaction_coords",
               "-r",
               radius,
               paste("-D fe_", format(radius, nsmall=6), sep=""),
               "-B nn",
               "-i network_end_node_traj.dat",
               paste("-o", output))

  message("Computing microstates from density network ...")
  run.cmds(cmds)

  if (sorted) {
    # reorder state names by descending population
    output <- paste(dir, output, sep="/")
    traj <- data.table::fread(output,
                              verbose=FALSE,
                              showProgress=FALSE)
    states <- unique(traj$V1)
    pops <- sapply(states, function(s){sum(traj$V1 == s)})

    # match old to new names
    matching <- data.frame(states, pops) %>%
                arrange(desc(pops))      %>%
                mutate(id=1:length(states))

    # convert to simple dictionary for fast retrieval
    dict <- list()
    for (s in states) {
      dict[s] <- matching %>%
        filter(states==s) %>%
        select(id)
    }
    # convert trajectory
    traj <- sapply(traj$V1, function(s){dict[[s]]})
    data.table::fwrite(list(traj),
                       output,
                       col.names=FALSE,
                       verbose=FALSE,
                       showProgress=FALSE)
  }
  message("... finished.")
}

#' Rename microstates.
#'
#' Rename microstates from 1 to N, sorted by descending populations.
#'
#' @param microstates Character or vector, either filename or microstate
#'  trajectory.
#' @param output Character, filename of renamed microstate trajectory.
#'  If \code{NULL} no output file is generated.
#' @return renamed microstate trajectory
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate count arrange
#' @export
clustering.microstates.renamed <- function(microstates, output=NULL) {
  if (is.character(microstates)) {
    microstates <- data.table::fread(microstates,
                                     verbose=FALSE,
                                     showProgress=FALSE)[[1]]
  }
  microstates_renamed <- microstates
  n_microstates <- length(unique(microstates))
  counts <- data.frame(state=microstates) %>%
            count(state)                  %>%
            arrange(desc(n))              %>%
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

#'
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

