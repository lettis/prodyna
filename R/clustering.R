#' Prepare the clustering output directory.
#'
#' The clustering output directory is created if it does not exist already.
#' A file named 'reaction_coords' that is linked to the reaction coordinates
#' file is created within that directoy.
#'
#' @param rc Character, name of the reaction coordinates file (normalized).
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#' @return absolute path to the output directory

clustering.prepareOutputDir <- function(rc, dir) {

  # default directory name
  if (is.null(dir)) {
    dir <- dirname(rc)
  }

  # create clustering directory
  if (!dir.exists(dir)) {
    message(paste("Creating output directory", dir))
    dir.create(dir)
  }

  normalizePath(dir)
}


#' Neighbourhood populations and free energies for density-based clustering.
#'
#' Compute neighbourhood populations and free energies for different radii.
#'
#' For each frame in the reaction coordinates file \code{rc} and each radius
#' in \code{radii} the neighbourhood population and corresponding free energies
#' are computed. The results are written to <prefix_pop>_<radius> and
#' <prefix_fe>_<radius>, respectively.
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radii Numeric vector, radii. The radius defines the size of the
#'  neighbourhood.
#' @param prefix_fe Character, free energy output file prefix
#'  (underscore and radius will be appended).
#' @param prefix_pop Character, populations output file prefix
#'  (underscore and radius will be appended).
#' @export
clustering.estimate.radii <- function(rc, radii, prefix_fe, prefix_pop) {

  if (!file.exists(rc)) {
    stop(msg("missingFile", rc))
  }

  args <- paste("density",
                "--file",
                rc,
                "--radii",
                paste(radii, collapse = " "),
                "--free-energy",
                prefix_fe,
                "--population",
                prefix_pop)

  message("Computing neighbourhood populations and free energies for density-based clustering ... ",
          appendLF=F)

  run.cmd(get.binary("clustering"), args)

  message("finished.")
}


#' Neighbourhood computation for efficient density-based geometric clustering.
#'
#' TODO: description
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radius Numeric, radius for density estimation.
#' @param prefix_fe Character, prefix of the free energy file (underscore and
#'  radius will be appended).
#' @param fname_nn Character, name of the nearest-neighbours output file.
#' @param prefixOnly Logical, is \code{fname_fe} specified in terms of a prefix
#'  only? (default: \code{TRUE}).
#' @export
clustering.compute.neighborhood <- function(rc, radius, prefix_fe, fname_nn) {

  # append radius to free energy file prefix
  fname_fe <- paste(prefix_fe, format(radius, nsmall=6), sep="_")

  # check if input files exist
  if(!file.exists(rc) || !file.exists(fname_fe)) {
    stop(msg("missingFile", c(rc, fname_fe)))
  }

  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                fname_fe,
                "--nearest-neighbors",
                fname_nn)

  message("Running neighborhood computation for density-based clustering ... ",
          appendLF=F)

  run.cmd(get.binary("clustering"), args)

  message("finished.")
}


#' Screening (density-based geometric clustering).
#'
#' Iterative clustering of frames below a specific free energy level.
#' The resulting state assignments are written to <prefix_clust>.<threshold>
#' where a state assignment of 0 indicates a frame with free energy larger
#' than the cut-off value.
#'
#' The parameters \code{min}, \code{max}, and \code{step} define which energy
#' thresholds are considered. The default value for both \code{min} and
#' \code{step} is 0.1.
#' If \code{max} is \code{NULL}, the maximum threshold corresponds to the
#' maximum free energy level in the dataset.
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radius Numeric, radius for density estimation.
#' @param prefix_fe Character, prefix of the free energy file (underscore and
#'  radius will be appended).
#' @param fname_nn Character, name of the nearest-neighbours file.
#' @param prefix_fe Character, prefix of the clustering output files (underscore and
#'  energy threshold will be appended).
#' @param min Numeric, minimum free energy threshold (default: 0.1)
#' @param max Numeric, maximum free energy threshold (default: \code{NULL},
#'  i.e. the maximum free energy level in the dataset)
#' @param step Numeric, step size by which the free energy threshold
#'  is increased.
#' @export
clustering.screening <- function(rc, radius, prefix_fe, fname_nn, prefix_clust, min=0.1, step=0.1, max=NULL) {

  feParams <- c(min, step, max)       # c(NULL, 1) == c(1)

  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                paste(prefix_fe, format(radius, nsmall=6), sep="_"),
                "--nearest-neighbors-input",
                fname_nn,
                "--threshold-screening",
                paste(feParams, collapse=" "),
                "--output",
                prefix_clust)

  message("Running screening ...")

  run.cmd(get.binary("clustering"), args)

  message("... finished.")
}


#' Construct density network from landscape screening.
#'
#' Output files are
#' \itemize{
#' \item remapped_clust_<feThreshold>
#' \item network_nodes.dat, network_links.datm, network_leaves.dat
#' \item network_end_node_traj.dat
#' \item network_visualization.html
#' }
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param minpop Numeric, minimum population to count microstate as node.
#' @param dir Character, name of the output directory. If \code{NULL} the
#'  output directory is assumed to be <rc>.clustering.
#'  The ouput directory is created if it does not exist already.
#' @param min Numeric, minimum free energy threshold
#' @param max Numeric, maximum free energy threshold
#' @param step Numeric, step size. In every iteration the free energy threshold
#'  is increased by this value.
#' @export
clustering.densityNetwork <- function(minpop, inputDir, outputDir=NULL, min=NULL, max=NULL, step=NULL) {

  ## TODO: min, max

  # dir <- clustering.prepareOutputDir(rc, dir)

  # construct density network
  cmds <- paste("cd", inputDir, ";")
  cmds <- paste(cmds,
                get.binary("clustering"),
                "network",
                "-p",
                minpop,
                "--step",
                step)

  message("Constructing density network ...")

  run.cmds(cmds)

  message("... finished.")
}


#' Compute Microstate trajectory.
#'
#' Assign frames to a distinct microstate using the tree structure obtained by
#' \code{\link{clustering.densityNetwork}}.
#' The resulting assignment is written to <output>.
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




