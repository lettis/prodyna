# List files in several directories
list.filesHelper <- function(dirs, p) {
  l <- as.list(unlist(lapply(dirs, function(d) {list.files(d, pattern = p, full.names = T)})))
  if (length(l) == 1 && length(l[[1]]) == 0)
    return(NULL)
  else
    return(l)
}

# List reference structure files (.pdb)
.list.reference <- function(dirs) {
  list.filesHelper(dirs, p = "*.pdb")
}

# List trajectory files (.xtc)
.list.traj <- function(dirs) {
  list.filesHelper(dirs, p = "*.xtc$")
}
# List files with dihedral angles (.dih)
.list.dihedrals <- function(dirs) {
  list.filesHelper(dirs, p = "*.dih$")
}

# List generic PCA files (cov, proj, etc.)
.list.PCA <- function(dirs, file_pattern) {
  PCA <- list()
  suffix <- c("proj", "cov", "vec", "val", "stats")
  suffix <- paste(rep(suffix, 2), c("", "n"), sep="")

  for (s in suffix) {
    PCA[[s]] <- list.filesHelper(dirs,
                                 p=paste(file_pattern, ".", s, "$", sep=""))
  }
  return(PCA)
}

# List specifically dPCA+ related files
.list.dPCAplus <- function(dirs) {
  .list.PCA(dirs, "*.dih")
}

# List C$_\alpha$ distance files
.list.caDists <- function(dirs) {
  list.files(dirs,
             pattern="*.dist_[[:digit:]]*_[[:digit:]]*$")
}

# List specifically C$_\alpha$-distance PCA related files
.list.caPCA <- function(dirs) {
  ca_dist_files <- .list.caDists(dirs)
  caPCA <- list()
  for (f in ca_dist_files) {
    suffix <- strsplit(f, ".", fixed=TRUE)[[1]]
    suffix <- suffix[length(suffix)]
    caPCA[[suffix]] <- .list.PCA(dirs, f)
  }
  return(caPCA)
}

# List files with filtered reaction coordinates
.list.reactionCoords <- function(dirs) {
  list.filesHelper(dirs, p = "reaction_coords")
}

#' Create a README file.
#'
#' Creates a README file to describe automatically generated files.
#' Note that the name of the reaction coordinates files needs to contain the
#' substring 'reaction_coords'.
#'
#' @param readme Character, name of the README file.
#' @param dirs Character, vector of directories that store the files to be
#'  listed in the README.
#' @importFrom rmarkdown render
#' @export
#'
make.readme <- function(readme  = "README", dirs) {

  sapply(dirs, normalizePath, USE.NAMES = F)

  sec <- function(title, level=1) {
    separator <- paste(rep("\n", 1), collapse="")
    indicator <- paste(rep("#", level), collapse="")
    paste(separator, paste(indicator, title), sep="")
  }
  txt <- function(text) {
    paste(text)
  }

  pd <- list()
  pd$ref <- .list.reference(dirs)
  pd$traj <- .list.traj(dirs)
  pd$dihedrals <- .list.dihedrals(dirs)
  pd$dPCAplus <- .list.dPCAplus(dirs)
  pd$caDists <- .list.caDists(dirs)
  pd$caPCA <- .list.caPCA(dirs)
  pd$reactionCoords <- .list.reactionCoords(dirs)

  proj_summary <- ""

  add_doc <- function(elem, doc) {
    if ((!is.null(elem)) && (!length(elem) == 0)) {
      proj_summary <<- c(proj_summary, doc)
    }
  }

  add_doc(pd$ref, c(sec("Reference Structure"),
                    txt(pd$ref)))
  add_doc(pd$traj, c(sec("Trajectory"),
                     txt(pd$traj)))
  add_doc(pd$dihedrals, c(sec("Dihedrals"),
                          txt(pd$dihedrals)))
  add_doc(pd$dPCAplus$proj, sec("dPCA+ Analysis"))
  add_doc(pd$dPCAplus$proj, c(sec("covariance-based", level=2),
                              txt(paste("**projection**",
                                        pd$dPCAplus$proj,
                                        "\n\n**covariance matrix**",
                                        pd$dPCAplus$cov,
                                        "\n\n**eigenvectors (columns)**",
                                        pd$dPCAplus$vec,
                                        "\n\n**eigenvalues**",
                                        pd$dPCAplus$val,
                                        "\n\n**stats**",
                                        pd$dPCAplus$stats))))
  add_doc(pd$dPCAplus$projn, c(sec("correlation-based", level=2),
                               txt(paste("**projection**",
                                         pd$dPCAplus$projn,
                                         "\n\n**correlation matrix**",
                                         pd$dPCAplus$covn,
                                         "\n\n**eigenvectors (columns)**",
                                         pd$dPCAplus$vecn,
                                         "\n\n**eigenvalues**",
                                         pd$dPCAplus$valn,
                                         "\n\n**stats**",
                                         pd$dPCAplus$statsn))))
  if (!is.null(pd$caDists) & length(pd$caDists) > 0) {
    proj_summary <- c(proj_summary, sec("C$_\\alpha$-distances"))
    for (d in pd$caDists) {
      add_doc(d, txt(d))
    }
  }
  if (!is.null(pd$caPCA) & length(pd$caPCA) > 0) {
    proj_summary <- c(proj_summary, sec("C$_\\alpha$-PCA"))
    for (i in 1:length(pd$caPCA)) {
      d <- pd$caPCA[i]
      proj_summary <- c(proj_summary, sec(names(d), level=2))
      d <- d[[1]]
      add_doc(d$proj, c(sec("covariance-based", level=3),
                        txt(paste("**projection**",
                                  d$proj,
                                  "\n\n**covariance matrix**",
                                  d$cov,
                                  "\n\n**eigenvectors (columns)**",
                                  d$vec,
                                  "\n\n**eigenvalues**",
                                  d$val,
                                  "\n\n**stats**",
                                  d$stats))))
      add_doc(d$projn, c(sec("correlation-based", level=3),
                         txt(paste("**projection**",
                                   d$projn,
                                   "\n\n**correlation matrix**",
                                   d$covn,
                                   "\n\n**eigenvectors (columns)**",
                                   d$vecn,
                                   "\n\n**eigenvalues**",
                                   d$valn,
                                   "\n\n**stats**",
                                   d$statsn))))
    }
  }
  if (length(pd$reactionCoords) > 0) {
    proj_summary <- c(proj_summary, sec("Reaction Coordinates"))
    lapply(pd$reactionCoords, function(fname) {
      if(substr(fname, nchar(fname)-5+1, nchar(fname)) != ".info" ) {
        proj_summary <<- c(proj_summary,
                           sec(fname, level="2"),
                           txt(
                             readLines(
                               paste(fname, "info", sep="."))))
      }
    })
  }

  cat(paste(c(proj_summary, ""), sep="", collapse="\n"),
      file=readme)
  opts <- list()
  opts["toc"] <- TRUE
  sink("/dev/null")
  render(readme, "pdf_document", output_options = opts, quiet = T)
  sink()
}

