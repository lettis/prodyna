

#' construct transition count matrix from trajectory
#'
#' Jumps are encoded as jump from row-index to column-index, i.e.
#' C[1,2] == 10 means 10 jumps from state 1 to state 2, while
#' C[2,1] == 20 means 20 jumps from state 2 to state 1.
#'
#' @param traj State trajectory, either encoded as vector or a filename.
#' @param lag time lag to use for counting jumps: traj[i] -> traj[i+lag].
#' @param traj_ids vector of same length as state trajectory. Assigns an id
#'                 to every frame identifying its trajectory, if state
#'                 trajectory is actually composed of several, independently run
#'                 simulations.
#' @export
msm.countMatrix <- function(traj, lag, traj_ids=NULL) {
  stopifnot(is.vector(traj) | is.character(traj))
  if (is.character(traj)) {
    traj <- data.table::fread(traj, verbose=FALSE)[[1]]
  }
  if (is.null(traj_ids)) {
    traj_ids <- rep(1, length(traj))
  }
  states <- sort(unique(traj))
  n_states <- length(states)
  C <- matrix(0, nrow=n_states, ncol=n_states)
  for (i in 1:(length(traj)-lag)) {
    if (traj_ids[i] == traj_ids[i+lag]) {
      s1 <- traj[i]
      s2 <- traj[i+lag]
      C[s1,s2] <- C[s1,s2] + 1
    }
  }
  C
}

#' construct transition matrix from trajectory
#' @param traj State trajectory, either encoded as vector or a filename.
#' @param lag time lag to use for counting jumps: traj[i] -> traj[i+lag].
#' @param traj_ids vector of same length as state trajectory. Assigns an id
#'                 to every frame identifying its trajectory, if state
#'                 trajectory is actually composed of several, independently run
#'                 simulations.
#' @param row.normalized if TRUE, transition matrix will be row-normalized, encoding jumps by from(row) to(col). Else (default) matrix is column-normalized encoding jumps from(col) to(row).
#' @param enforce.db Enforce detailed balance by averaging over
#'                   forward and backward rates.
#' @export
msm.transitionMatrix <- function(traj,
                                 lag,
                                 traj_ids=NULL,
                                 row.normalized=FALSE) {
  C <- msm.countMatrix(traj, lag, traj_ids)
  for (i in 1:nrow(C)) {
    C[i,] <- C[i,] / sum(C[i,])
  }
  if ( ! row.normalized) {
    C <- t(C)
  }
  C
}


#' compute transition- and life-times from transition matrix
#' @param tmat transition matrix.
#' @param lag lagtime (expressed in number of frames) used for tmat estimation.
#' @param dt time step for a single frame. defines unit of time for output.
#' @param row.normalized Is tmat row-normalized, i.e. transitions are encoded
#'                       as  T_ij: transition i -> j (default: FALSE) or
#'                       col-normalized (T_ij: j -> i, this is the default)?
#' @export
msm.transitionTimes <- function(tmat,
                                lag,
                                dt,
                                row.normalized=FALSE) {
  if (row.normalized) {
    tmat <- t(tmat)
  }
  # transition times
  ttimes <- (dt*lag) / tmat
  N <- nrow(tmat)
  for (i in 1:N) {
    # life times
    ttimes[i,i] <- 1 / sum(tmat[(1:N)!=i, i] / (dt*lag))
  }
  if (row.normalized) {
    ttimes <- t(ttimes)
  }

  ttimes
}




#' simulate state trajectory from Markov State Model
#'
#' Generates a state trajectory of length n_steps+1 (n_steps simulation steps plus
#' the initial state). State names are in the range of 1:nrow(T).
#'
#' @param T column-normalized transition matrix encoding the MSM.
#' @param initial_state State number of state to start with. Must be in range 1:nrow(T).
#' @param n_steps Number of simulation steps.
#' @export
msm.sim <- function(T, initial_state, n_steps){
  n_states <- nrow(T)
  traj <- rep(initial_state, n_steps)
  for (i in 2:n_steps) {
    traj[i] <- sample.int(n_states,
                          size=1,
                          prob=T[1:n_states,traj[i-1]])
  }
  traj
}

#' compute mean first passage times (MFPT)
#'
#' Compute MFPT from one state to all other states.
#' Computation is done by Markov Chain Monte Carlo sampling of the
#' given (col-normalized) transition matrix.
#'
#' @param T column-normalized transition matrix encoding the MSM.
#' @param initial_state State number of state for which MFPTs should be
#'        determined.
#' @param n_chains Number of chains to be run (more chains = better sampling).
#' @param n_steps Number of timesteps. Should be higher then the slowest
#'        timescale of the system in order to sample of MFPTs.
#' @export
msm.mfpt <- function(T, initial_state, n_chains, n_steps) {
  n_states <- nrow(T)
  mfpt <- matrix(NA, nrow=n_chains, ncol=n_states)
  for (i in 1:n_chains) {
    traj <- msm.sim(T, initial_state, n_steps)
    mfpt[i,] <- sapply(1:n_states, function(j) {
      match(j, traj)
    })
  }
  sapply(1:n_states, function(j) {
    mean(mfpt[,j], na.rm=T)
  })
}
