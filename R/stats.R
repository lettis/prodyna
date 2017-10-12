
#' Autocorrelation.
#'
#' Correlation between a signal and itself delayed in time.
#'
#' @param x Character or data.frame, either a filename or a data frame
#'  providing the data.
#' @param lag.max Numeric, maximal correlation lag.
#'  If < 1, expressed in fractions of total length (default: 0.25).
#' @param type Character, one of \code{c("correlation", "covariance")}. NOT USED
#' @param columns Numeric vector, selected observables for which the
#'  autocorrelation should be computed. If \code{NULL} compute all.
#' @param circular Logical, should autocorrelation be computed for circular
#'  variables, e.g. dihedral angles? (default: \code{FALSE})
#' @param circ.units Character, circular units, one of \code{c("deg", "rad")}
#' @importFrom data.table fread
#' @importFrom NISTunits NISTdegTOradian
#' @return data.frame where row i corresponds to the autocorrelation for a delay
#'  of i-1 frames.
#' @export
stats.autocor <- function(x, lag.max=0.25, type="correlation", columns=NULL, circular=FALSE, circ.units="deg") {
  # prepare data import
  if (is.character(x)) {                    # interpret as filename
    colSelect <- function(col) {
      fread(input=x, select=col, verbose=F, showProgress=F)$V1
    }
    if (is.null(columns)) {
      columns <- names(fread(x, nrows=1))
    }
  } else {                                  # interpret as data.frame
    colSelect <- function(col) { x[[col]] }
    if (is.null(columns)) {
      columns <- 1:ncol(x)
    }
  }
  # compute the autocorrelations
  cors <- lapply(columns, function(c) {
      xs <- colSelect(c)
      N  <- length(xs)
      if (lag.max < 1) {
        tau_max <- round(lag.max*N)
      } else {
        tau_max <- lag.max
      }
      if (circular) {
        # ensure radians
        if (tolower(circ.units) == "deg") {
          xs <- NISTdegTOradian(xs)
        }
        # compute ACF (cov) on circular data
        rho <- circacf(xs, tau_max)
      } else {
        # compute ordinary ACF, but with parallelized code
        rho <- fastacf(xs, tau_max)
      }
      # normalize -> ACF (corr)
      rho/rho[[1]]
  })

  names(cors) <- columns
  return(data.frame(cors))
}
