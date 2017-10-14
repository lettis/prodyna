

#' Transition count matrix.
#'
#' Construct a transition count matrix from a trajectory.
#' Jumps are encoded as jump from row-index to column-index, i.e.
#' C[1,2] == 10 means 10 jumps from state 1 to state 2, while
#' C[2,1] == 20 means 20 jumps from state 2 to state 1.
#'
#' @param traj Character or numeric vector, name of trajectory file or vector
#'  describing the state trajectory.
#' @param lag Numeric, time lag to use for counting jumps expressed in number
#'  number of frames, traj[i] -> traj[i+lag].
#' @param traj_ids Numeric vector, frame ids identifying the trajectory. These
#'  should be set if the state trajectory is actually composed of several,
#'  independently run simulations. If \code{NULL} a single trajectory is assumed.
#'  Must be of same length as the state trajectory.
#' @export
msm.countMatrix <- function(traj, lag, traj_ids=NULL) {
  stopifnot(is.vector(traj) | is.character(traj))
  if (is.character(traj)) {
    traj <- data.table::fread(traj, verbose=FALSE, showProgress=FALSE)[[1]]
  }
  if (is.null(traj_ids)) {
    traj_ids <- rep(1, length(traj))
  }

  states   <- sort(unique(traj))
  n_states <- length(states)

  C <- matrix(0, nrow=n_states, ncol=n_states)

  for (i in 1:(length(traj)-lag)) {
    if (traj_ids[i] == traj_ids[i+lag]) {
      s1 <- traj[i]
      s2 <- traj[i+lag]
      C[s1,s2] <- C[s1,s2] + 1
    }
  }
  return(C)
}

#' Transition Matrix.
#'
#' Construct a transition matrix from a trajectory.
#'
#' @param traj Character or numeric vector, name of trajectory file or vector
#'  describing the state trajectory.
#' @param lag Numeric, time lag to use for counting jumps expressed in number
#'  number of frames, traj[i] -> traj[i+lag].
#' @param traj_ids Numeric vector, frame ids identifying the trajectory. These
#'  should be set if the state trajectory is actually composed of several,
#'  independently run simulations. If \code{NULL} a single trajectory is assumed.
#'  Must be of same length as the state trajectory.
#' @param row.normalized Logical, if \code{FALSE} (default) the transition
#'  matrix is column-normalized, i.e. the entry T_ij encodes the probability
#'  of a jump from j to i.
#'  Otherwise, transition matrix will be row-normalized, i.e. the entry T_ij
#'  encodes the probability of a jump from i to j.
#' @export
msm.transitionMatrix <- function(traj, lag, traj_ids=NULL, row.normalized=FALSE) {
  C <- msm.countMatrix(traj, lag, traj_ids)

  for (i in 1:nrow(C)) {
    C[i,] <- C[i,] / sum(C[i,])
  }
  if (!row.normalized) {
    return(t(C))
  } else {
    return(C)
  }
}


#' Transition times and life times.
#'
#' Compute transition and life times from a transition matrix.
#'
#' @param tmat Matrix, transition matrix.
#' @param lag Numeric, lagtime (expressed in number of frames) that was used
#'  for estimating the transition matrix estimation.
#' @param dt Numeric, time step for a single frame. This defines the unit of
#'  time for the output.
#' @param row.normalized Logical, Is tmat row-normalized, i.e. transitions are encoded
#'                       as  T_ij: transition i -> j (default: FALSE) or
#'                       col-normalized (T_ij: j -> i, this is the default)?
#' @param row.normalized Logical, if \code{FALSE} (default) the transition
#'  matrix \code{tmat} is column-normalized, i.e. the entry \code{tmat[i,j]}
#'  encodes the probability of a jump from j to i.
#'  Otherwise, transition matrix will be row-normalized, i.e. the entry
#'  \code{tmat[i,j]} encodes the probability of a jump from i to j.
#' @export
msm.transitionTimes <- function(tmat, lag, dt, row.normalized=FALSE) {
  if (row.normalized) {
    tmat <- t(tmat)
  }
  # transition times
  ttimes <- (dt*lag)/tmat
  N      <- nrow(tmat)

  for (i in 1:N) {
    # life times
    ttimes[i,i] <- 1/sum(tmat[(1:N)!=i, i]/(dt*lag))
  }
  if (row.normalized) {
    return(t(ttimes))
  } else {
    return(ttimes)
  }
}


#' MSM State Trajectory Simulation.
#'
#' Generate a state trajectory from a Markov State Model. The state trajectory
#' will have length \code{n_steps+1} (n_steps simulation steps plus the
#' initial state).
#' State names are in the range of \code{1:nrow(T)}.
#'
#' @param T Matrix, column-normalized transition matrix encoding the MSM.
#' @param initial_state Numeric, initial state (must be in range
#'  \code{1:nrow(T)}).
#' @param n_steps Numeric, number of simulation steps.
#' @export
msm.sim <- function(T, initial_state, n_steps){
  n_states <- nrow(T)
  traj <- rep(initial_state, n_steps)
  for (i in 2:n_steps) {
    traj[i] <- sample.int(n_states,
                          size=1,
                          prob=T[1:n_states,traj[i-1]])
  }
  return(traj)
}

#' Mean first passage times (MFPT)
#'
#' Compute the Mean first passage times (MFPT) from one state to all other
#' states.
#' Computation is done by Markov Chain Monte Carlo sampling of the
#' given (col-normalized) transition matrix.
#'
#' @param T Matrix, column-normalized transition matrix encoding the MSM.
#' @param initial_state Numeric, state for which MFPTs should be determined.
#' @param n_chains Numeric, number of chains to be run
#'  (more chains = better sampling).
#' @param n_steps Numeric, number of timesteps. Should be higher then the
#'  slowest timescale of the system in order to sample of MFPTs.
#' @export
msm.mfpt <- function(T, initial_state, n_chains, n_steps) {
  n_states <- nrow(T)
  mfpt     <- matrix(NA, nrow=n_chains, ncol=n_states)

  for (i in 1:n_chains) {
    traj     <- msm.sim(T, initial_state, n_steps)
    mfpt[i,] <- sapply(1:n_states, function(j) {match(j, traj)})
  }
  colMeans(mfpt, na.rm=TRUE)
}
