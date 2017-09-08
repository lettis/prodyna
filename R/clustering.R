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
#' @param prefix_clust Character, prefix of the clustering output files (underscore and
#'  energy threshold will be appended).
#' @param min Numeric, minimum free energy threshold (default: 0.1)
#' @param max Numeric, maximum free energy threshold (default: \code{NULL},
#'  i.e. the maximum free energy level in the dataset)
#' @param step Numeric, step size by which the free energy threshold
#'  is increased.
#' @export
clustering.screening <- function(rc, radius, prefix_fe, fname_nn, prefix_clust, min=0.1, max=NULL, step=0.1) {

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

  message("Running screening ... ", appendLF=F)

  run.cmd(get.binary("clustering"), args)

  message("finished.")
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
#' They are written to the directory of the <prefix_clust> files.
#'
#' TODO clustering files currently need to be named clust.<radius>
#'
#' @param minpop Numeric, minimum population to count a microstate as node.
#' @param prefix_clust Character, prefix of the clustering files.
#' @param step Numeric, step size, i.e. free energy differences between nodes.
#' @export
clustering.densityNetwork <- function(minpop, prefix_clust, step=0.1) {

  dir <- dirname(prefix_clust)
  prefix_clust <- basename(prefix_clust)

  cmds <- paste("cd", dir, ";",
                get.binary("clustering"),
                "network",
                #"--basename",
                #prefix_clust,
                "--minpop",
                minpop,
                "--step",
                step)

  message("Constructing density network ... ", appendLF=F)

  run.cmds(cmds)

  message("finished.")
}


#' Compute Microstate trajectory.
#'
#' Assign frames to a distinct microstate using the tree structure obtained by
#' \code{\link{clustering.densityNetwork}}.
#' The resulting assignment is written to file <microstates>.
#'
#' @param rc Character, name of reaction coordinates file.
#' @param radius Numeric, density estimation radius.
#' @param prefix_fe Character, prefix of the free energy files (underscore and
#'  radius will be appended).
#' @param nn Character, name of the nearest-neighbours file.
#' @param init Character, name of the file specifying the initial state
#'  definition.
#' @param microstates Character, name of the file storing the resulting state
#'  trajectory.
#' @param sorted Logical, should microstate names be ordered by descending
#'  populations? (Default: \code{TRUE})
#' @importFrom dplyr desc arrange mutate filter select
#' @importFrom magrittr "%>%"
#' @export
clustering.microstates <- function(rc, radius, prefix_fe, nn, init, microstates, sorted=TRUE) {

  # run neighborhood computation
  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                paste(prefix_fe, format(radius, nsmall=6), sep="_"),
                "--nearest-neighbors-input",
                nn,
                "--input",
                init,
                "--output",
                microstates)

  message("Computing microstates from density network ... ", appendLF=F)
  run.cmd(get.binary("clustering"), args)

  if (sorted) {
    # reorder state names by descending population
    traj   <- data.table::fread(microstates, verbose=FALSE, showProgress=FALSE)
    states <- unique(traj$V1)
    pops   <- sapply(states, function(s){sum(traj$V1 == s)})

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
                       microstates,
                       col.names=FALSE,
                       verbose=FALSE,
                       showProgress=FALSE)

  }
  message("finished.")
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




