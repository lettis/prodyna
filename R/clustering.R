
#' Neighbourhood populations and free energies for density-based clustering.
#'
#' Compute neighbourhood populations and free energies for different radii.
#'
#' For each frame in the reaction coordinates file \code{rc} and each radius
#' in \code{radii} the neighbourhood population and corresponding free energies
#' are computed. The results are written to <pop_prefix>_<radius> and
#' <fe_prefix>_<radius>, respectively.
#'
#' @param rc Character, name of the reaction coordinates file.
#' @param radii Numeric vector, radii. The radius defines the size of the
#'  neighbourhood.
#' @param fe_prefix Character, free energy output file prefix
#'  (underscore and radius will be appended). Default is "fe".
#' @param pop_prefix Character, populations output file prefix
#'  (underscore and radius will be appended). Default is "pop".
#' @export
clustering.estimate.radii <- function(rc, radii, fe_prefix, pop_prefix) {

  if (!file.exists(rc)) {
    stop(msg("missingFile", rc))
  }

  args <- paste("density",
                "--file",
                rc,
                "--radii",
                paste(radii, collapse = " "),
                "--free-energy",
                fe_prefix,
                "--population",
                pop_prefix)

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
#' @param fe_prefix Character, prefix of the free energy file (underscore and
#'  radius will be appended). Default is "fe".
#' @param nn Character, name of the nearest-neighbours output file.
#' @export
clustering.compute.neighborhood <- function(rc, radius, fe_prefix, nn) {

  # append radius to free energy file prefix
  fe_fname <- paste(fe_prefix, format(radius, nsmall=6), sep="_")

  # check if input files exist
  if(!all(file.exists(c(rc, fe_fname)))) {
    stop(msg("missingFile", c(rc, fe_fname)))
  }

  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                fe_fname,
                "--nearest-neighbors",
                nn)

  message("Running neighborhood computation for density-based clustering ... ",
          appendLF=F)

  run.cmd(get.binary("clustering"), args)

  message("finished.")
}


#' Screening (density-based geometric clustering).
#'
#' Iterative clustering of frames below a specific free energy level.
#' The resulting state assignments are written to <clust_prefix>.<threshold>
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
#' @param fe_prefix Character, prefix of the free energy file (underscore and
#'  radius will be appended).
#' @param nn Character, name of the nearest-neighbours file.
#' @param clust_prefix Character, prefix of the clustering output files
#'  (underscore and energy threshold will be appended).
#' @param min Numeric, minimum free energy threshold (default: 0.1)
#' @param max Numeric, maximum free energy threshold (default: \code{NULL},
#'  i.e. the maximum free energy level in the dataset)
#' @param step Numeric, step size by which the free energy threshold
#'  is increased.
#' @export
clustering.screening <- function(rc, radius, fe_prefix, nn, clust_prefix, min=0.1, max=NULL, step=0.1) {

  feParams <- c(min, step, max)       # c(NULL, 1) == c(1)

  # append radius to free energy file prefix
  fe_fname <- paste(fe_prefix, format(radius, nsmall=6), sep="_")

  # check if input files exist
  if(!all(file.exists(rc, fe_fname, nn))) {
    stop(msg("missingFile", c(rc, fe_fname, nn)))
  }

  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                fe_fname,
                "--nearest-neighbors-input",
                nn,
                "--threshold-screening",
                paste(feParams, collapse=" "),
                "--output",
                clust_prefix)

  message("Running screening ... ", appendLF=F)

  run.cmd(get.binary("clustering"), args)

  message("finished.")
}


#' Construct density network from landscape screening.
#'
#' Output files are
#' \itemize{
#' \item remapped_clust_<feThreshold>
#' \item network_nodes.dat, network_links.dat, network_leaves.dat
#' \item network_end_node_traj.dat
#' \item network_visualization.html
#' }
#' Note that these files will be written to the directory that contains the
#' <clust_prefix>_<radius> files.
#'
#' TODO clustering files currently need to be named clust.<radius>
#'
#' @param minpop Numeric, minimum population to count a microstate as node.
#' @param clust_prefix Character, prefix of the clustering files obtained by
#'  screening the free energy landscape (se \code\link{clustering.screening}).
#' @param step Numeric, step size, i.e. free energy differences between nodes.
#' @export
clustering.densityNetwork <- function(minpop, clust_prefix="clust", step=0.1) {

  dir          <- normalizePath(dirname(clust_prefix))
  clust_prefix <- basename(clust_prefix)

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
#' \link{clustering.densityNetwork}.
#' The resulting assignment is written to file <microstates>.
#'
#' @param rc Character, name of reaction coordinates file.
#' @param radius Numeric, density estimation radius.
#' @param fe_prefix Character, prefix of the free energy files (underscore and
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
#' @importFrom data.table fwrite fread
#' @export
clustering.microstates <- function(rc, radius, fe_prefix, nn, init, microstates, sorted=TRUE) {

  # append radius to free energy file prefix
  fe_fname <- paste(fe_prefix, format(radius, nsmall=6), sep="_")

  # check if input files exist
  if(!all(file.exists(rc, fe_fname, nn, init))) {
    stop(msg("missingFile", c(rc, fe_fname, nn, init)))
  }

  # run neighborhood computation
  args <- paste("density",
                "--file",
                rc,
                "--radius",
                radius,
                "--free-energy-input",
                fe_fname,
                "--nearest-neighbors-input",
                nn,
                "--input",
                init,
                "--output",
                microstates)

  message("Computing microstates from density network ... ", appendLF=F)
  run.cmd(get.binary("clustering"), args)
  message("finished.")

  if (sorted) {
    # reorder state names by descending population
    traj   <- fread(microstates, verbose=FALSE, showProgress=FALSE)
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

    fwrite(list(traj), microstates, col.names=F, verbose=F, showProgress=F)
  }

}

#' Rename microstates.
#'
#' Rename microstates from 1 to N, sorted by descending populations.
#'
#' @param microstates Character or vector, either filename or microstate
#'  trajectory.
#' @param output Character, filename of renamed microstate trajectory.
#'  If \code{NULL}, no output file is created.
#' @return renamed microstate trajectory
#' @importFrom magrittr "%>%"
#' @importFrom data.table fwrite fread
#' @importFrom dplyr mutate count arrange
#' @export
clustering.microstates.renamed <- function(microstates, output=NULL) {
  if (is.character(microstates)) {
    microstates <- fread(microstates, verbose=FALSE, showProgress=FALSE)[[1]]
  }
  microstates_renamed <- microstates
  n_microstates       <- length(unique(microstates))
  counts <- data.frame(state=microstates) %>%
            count(state)                  %>%
            arrange(desc(n))              %>%
            mutate(idx=1:n_microstates)

  for (state in unique(microstates)) {
    microstates_renamed[microstates==state] <- counts[counts$state==state,]$idx
  }
  if (!is.null(output)) {
    fwrite(data.frame(microstates_renamed),
           sep=" ",
           col.names    = FALSE,
           verbose      = FALSE,
           showProgress = FALSE)
  }
  return(microstates_renamed)
}


#' MPP dynamic clustering
#'
#' Peform dynamic clustering using the most-probable-path algorithm to assign
#' geometric microstates to metastable states.
#' Either specify a single metastability criterion \code{Q} or a range of
#' metatstability values by specifying \code{Qmin}, \code{Qmax}, \code{Qstep}.
#' (If \code{Q} is not \code{NULL}, \code{Qmin}, \code{Qmax}, \code{Qstep} will
#' be ignored.)
#'
#' @param microstates Character, name of the file storing the state trajectory
#'  obtained via geometric clustering.
#' @param fe_prefix Character, prefix of the free energy files (underscore and
#'  radius will be appended).
#' @param radius Numeric, density estimation radius.
#' @param lagtime Numeric, number of frames to be skipped when calculating.
#'  transition probabilites (timescales below this lagtime will be discarded).
#' @param mpp Character, prefix of the MPP population and MPP trajectory files.
#' @param Q Numeric, metastability criterion (\eqn{\in (0,1]})
#' @param Qmin Numeric, minimum metastability (\eqn{\in (0,1]})
#' @param Qstep Numeric, step size
#' @param Qmax Numeric, maximum metastability (\eqn{\in (0,1]})
#' @examples
#'  clustering.mpp(microstates = "microstates",
#'                 fe_prefix   = "freeEnergy",
#'                 radius      = 0.3,
#'                 lagtime     = 50,
#'                 mpp_prefix  = "mpp",
#'                 Q           = 0.3)
#'  clustering.mpp(microstates = "microstates",
#'                 fe_prefix   = "freeEnergy",
#'                 radius      = 0.3,
#'                 lagtime     = 50,
#'                 mpp_prefix  = "mpp",
#'                 Qmin        = 0.1,
#'                 Qstep       = 0.01,
#'                 Qmax        = 1)
#' @export
clustering.mpp <- function(microstates, fe_prefix, radius, lagtime, mpp_prefix, Q=NULL, Qmin=0.01, Qstep=0.01, Qmax=1) {

  if (!all(c(Q, Qmin, Qstep, Qmax) > 0) || !all(c(Q, Qmin, Qstep, Qmax) <= 1)) {
    stop("Q, Qmin, Qstep, Qmax must all be within (0,1].")
  }
  # append radius to free energy file prefix
  fe_fname <- paste(fe_prefix, format(radius, nsmall=6), sep="_")

  if(!all(file.exists(microstates, fe_fname))) {
    stop(msg("missingFile", c(microstates, fe_fname)))
  }

  if(!is.null(Q)) {
    Qmin <- Q
    Qmax <- Q
  }
  # run mpp clustering
  args <- paste("mpp",
                "--input",
                microstates,
                "--free-energy-input",
                fe_fname,
                "--lagtime",
                lagtime,
                "--basename",
                mpp_prefix,
                "--qmin-from",
                Qmin,
                "--qmin-to",
                Qmax,
                "--qmin-step",
                Qstep
                )

  message("Running dynamic clutering using the MPP method ... ", appendLF=F)
  run.cmd(get.binary("clustering"), args)
  message("finished.")
}

#' Waiting Time Distributions
#'
#'  Compute waiting time distributions for the give window sizes.
#'  The resulting distributions will be written to
#'  <wtDistribution_prefix>_<window size>_<state id>.
#'
#' @param traj Character, name of the file holding the state trajectory.
#' @param wtd_prefix Character, prefix of the output file holding
#'  the waiting time distribution (window size and state id will be appended).
#' @param wsizes Numeric, vector of window size.
#' @export
clustering.wtDistributions <- function(traj, wtd_prefix, wsizes) {

  if(!file.exists(traj)) {
    stop(msg("missingFile", traj))
  }

  win_fname <- paste(traj, "windows", sep="_")
  unlink(win_fname)

  for(w in wsizes) {

    write(paste("*", w), file = win_fname)

    wtDistribution_fname <- paste(wtd_prefix, w, sep="_")

    args <- paste("coring",
                  "--states",
                  traj,
                  "--windows",
                  win_fname,
                  "--distribution",
                  wtDistribution_fname)

    run.cmd(get.binary("clustering"), args)

    unlink(win_fname)
  }
}

#' Variable Dynamic Coring
#'
#' A state transition is accepted only if the trajectory spends a minimum
#' time (given by the window sizes) in the new states. Otherwise the frame is
#' reassigned to the previous state. The resulting trajectory exhibits clear-cut
#' state boundaries.
#'
#' @param traj Character, name of the file holding the state trajectory.
#' @param traj_cored Character, name of the output file holding the cored
#'  state trajectory.
#' @param states Numeric, vector of state ids.
#' @param wsizes Numeric, vector of window sizes per state, i.e. use a coring
#'  window of size \code{wsizes[i]} for state \code{states[i]}.
#' @param wsize_default Numeric, window size that should be used for all states
#'  not in \code{states}.
#' @export
clustering.coring <- function(traj, traj_cored, states, wsizes, wsize_default=NULL) {

  if(!file.exists(traj)) {
    stop(msg("missingFile", traj))
  }

  win_fname <- paste(traj, "windows", sep="_")
  unlink(win_fname)

  if(!is.null(wsize_default)) {
    w_specification <- paste("*", wsize_default, "\n")
  } else {
    w_specification <- ""
  }

  for(i in length(wsizes)) {

    w_specification <- paste(w_specification,
                             states[i],
                             wsizes[i],
                             "\n")
  }

  write(w_specification, file=win_fname)

  args <- paste("coring",
                "--states",
                traj,
                "--windows",
                win_fname,
                "--output",
                traj_cored)

  run.cmd(get.binary("clustering"), args)

  unlink(win_fname)
}


