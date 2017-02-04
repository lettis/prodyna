

#' compute autocorrelations
#' @param x The data. Either a filename or a data frame with values.
#' @param lag.max The maximal correlation lag. If < 1, expressed in fractions of total length.
#' @param type One of c("correlation", "covariance")
#' @param columns The selected observables for which the autocorrelation should be computed. If NULL: compute all.
#' @param circular Compute autocorrelation for circular variables (e.g. dihedral angles).
#' @param circ.units Circular units, one of c("deg", "rad")
#' @export
stats.autocor <- function(x, lag.max=0.25, type="correlation", columns=NULL, circular=FALSE, circ.units="deg") {
  # prepare data import
  if (is.character(x)) {
    # interpret as filename
    colSelect <- function(col) {
      data.table::fread(input=x, select=col, verbose=FALSE, showProgress=FALSE)
    }
    if (is.null(columns)) {
      columns <- names(data.table::fread(x, nrows=1))
    }
  } else {
    # interpret as object
    colSelect <- function(col) {
      x$col
    }
    if (is.null(columns)) {
      columns <- names(x)
    }
  }
  # compute the autocorrelations
  cors <- lapply(columns, function(c) {
    xs <- colSelect(c)[[1]]
    N <- length(xs)
    if (lag.max < 1) {
      tau_max <- round(lag.max*N)
    } else {
      tau_max <- lag.max
    }
    if (circular) {
      # ensure radians
      if (tolower(circ.units) == "deg") {
        xs <- NISTunits::NISTdegTOradian(xs)
      }
      # compute ACF (cov) on circular data
      rho <- circacf(xs, tau_max)
    } else {
      # compute ordinary ACF, but with parallelized code
      rho <- fastacf(xs, tau_max)
    }
    # normalize -> ACF (corr)
    rho / rho[[1]]
  })

  cors
}
