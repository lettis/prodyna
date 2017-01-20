

#' construct transition count matrix from trajectory
#'
#' Jumps are encoded as jump from row-index to column-index, i.e.
#' C[1,2] == 10 means 10 jumps from state 1 to state 2, while
#' C[2,1] == 20 means 20 jumps from state 2 to state 1.
#'
#' @param traj State trajectory, either encoded as vector or a filename.
#' @param lag time lag to use for counting jumps: traj[i] -> traj[i+lag].
#' @export
msm.countMatrix <- function(traj, lag) {
  stopifnot(is.vector(traj) | is.character(traj))
  if (is.character(traj)) {
    traj <- data.table::fread(traj, verbose=FALSE)[[1]]
  }
  states <- sort(unique(traj))
  n_states <- length(states)
  C <- matrix(0, nrow=n_states, ncol=n_states)
  for (i in 1:(length(traj)-lag)) {
    s1 <- traj[i]
    s2 <- traj[i+lag]
    C[s1,s2] <- C[s1,s2] + 1
  }
  C
}

#' construct transtion matrix from trajectory
#' @param traj State trajectory, either encoded as vector or a filename.
#' @param lag time lag to use for counting jumps: traj[i] -> traj[i+lag].
#' @param row.normalized if TRUE, transition matrix will be row-normalized, encoding jumps by from(row) to(col). Else (default) matrix is column-normalized encoding jumps from(col) to(row).
#' @export
msm.transitionMatrix <- function(traj, lag, row.normalized=FALSE) {
  C <- msm.countMatrix(traj, lag)
  for (i in 1:nrow(C)) {
    C[i,] <- C[i,] / sum(C[i,])
  }
  if ( ! row.normalized) {
    C <- t(C)
  }
  C
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

  propagate <- function(state) {
    candidate <- sample(1:n_states, 1)
    while(runif(1) > T[candidate,state]) {
      candidate <- sample(1:n_states, 1)
    }
    candidate
  }

  traj <- initial_state
  for (i in 1:n_steps) {
    traj <- c(traj, propagate(traj[i]))
  }
  traj
}


