
#' read dihedrals
#' @return data.frame with \eqn{\phi} and \eqn{\psi} angles
#' @param resno Residue number. If NULL: read dihedrals from all residues.
#' @export
read.dihedrals <- function(resno=NULL) {
  .check.projectPath()
  # get project information
  pd <- projectInfo()
  if (is.null(resno)) {
    dih <- data.table::fread(get.fullPath(pd$dihedrals),
                             verbose=FALSE,
                             showProgress=FALSE)
    n_dih <- dim(dih)[2]/2
    res_ndx <- do.call("c", lapply(2:(n_dih+1), function(i){rep(i,2)}))
  } else {
    # TODO resno as vector?
    dih <- data.table::fread(get.fullPath(pd$dihedrals),
                             select=c(2*(resno-1)-1, 2*(resno-1)),
                             verbose=FALSE,
                             showProgress=FALSE)
    res_ndx <- c(resno, resno)
  }
  colnames(dih) <- paste(c("phi", "psi"), res_ndx, sep="")

  dih
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
