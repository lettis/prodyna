
#' Read dihedrals.
#'
#' Read the columns corresponding to the \eqn{\phi} and \eqn{\psi} angles of the
#' given residues.
#' Assumes first two columns correpond to dihedral angles of residue 2.
#'
#' @param resnos Numeric vector, residue numbers for which dihedral angles
#'  should be read. If \code{NULL} (default) read dihedrals for all residues.
#' @param dihedrals Character, name of the .dih file containing the dihedral
#'  angles.
#' @return A data.frame with \eqn{\phi} and \eqn{\psi} angles.
#' @export
read.dihedrals <- function(dihedrals, resnos=NULL, reslabel=NULL) {
  if (is.null(resnos)) {
    dih <- data.table::fread(dihedrals,
                             verbose=FALSE,
                             showProgress=FALSE)
    n_dih <- dim(dih)[2]/2
    resnos <- 2:(n_dih+1)
  } else {
    cols <- do.call("c", lapply(resnos, function(i){c(2*(i-1)-1, 2*(i-1))}))
    dih <- data.table::fread(dihedrals,
                             select=cols,
                             verbose=FALSE,
                             showProgress=FALSE)
  }
  res_ndx <- do.call("c", lapply(resnos, function(i){rep(i,2)}))
  colnames(dih) <- paste(c("phi", "psi"), res_ndx, sep="")

  return(dih)
}

#' cumulative fluctuations for given PCA
#' @param pca Selected PCA-procedure
#' @param corr Use correlation-based PCA
#' @return cumulative fluctuations
#' @examples
#' p <- prodyna::projectInfo()
#' prodyna::read.cumFlucts(p$dPCAplus)
#' prodyna::read.cumFlucts(p$caPCA$dist_3_6, corr=TRUE)
#' @export
read.cumFlucts <- function(pca, corr=FALSE) {
  .check.projectPath()
  #TODO test if exists, else run
  if (corr) {
    vals <- data.table::fread(pca$valn)
  } else {
    vals <- data.table::fread(pca$val)
  }
  vals <- vals$V1

  cumsum(vals/max(vals))
}
