
#' Project cache
#'
#' Holds all project-related information.
#' Automatically updated by generator functions.
.project_cache <- new.env()

#' get project information
#'
#' @export
project <- function() {
  if (exists("project",
             envir=.project_cache,
             inherits=FALSE)) {
    get("project", envir=.project_cache)
  } else {
    warning(.warnings$no_init)
    NULL
  }
}

## project file detection

#' List reference structure files (.pdb)
.list.reference <- function() {
  list.files(pattern = "*.pdb")
}
#' List trajectory files (.xtc)
.list.traj <- function() {
  list.files(pattern = "*.xtc$")
}
#' List files with dihedral angles (.dih)
.list.dihedrals <- function() {
  list.files(pattern = "*.dih$")
}
#' List generic PCA files (cov, proj, etc.)
.list.PCA <- function(file_pattern) {
  PCA <- list()
  suffix <- c("proj", "cov", "vec", "val", "stats")
  suffix <- paste(rep(suffix, 2), c("", "n"), sep="")
  for (s in suffix) {
    PCA[[s]] <- list.files(pattern=paste(file_pattern,
                                         ".",
                                         s,
                                         "$",
                                         sep=""))
  }

  PCA
}
#' List specifically dPCA+ related files
.list.dPCAplus <- function() {
  .list.PCA("*.dih")
}
#' List C$_\alpha$ distance files
.list.caDists <- function() {
  list.files(pattern="*.dist_[[:digit:]]*_[[:digit:]]*$")
}
#' List specifically C$_\alpha$-distance PCA related files
.list.caPCA <- function() {
  ca_dist_files <- .list.caDists()
  caPCA <- list()
  for (f in ca_dist_files) {
    suffix <- strsplit(f, ".", fixed=TRUE)[[1]]
    suffix <- suffix[length(suffix)]
    caPCA[[suffix]] <- .list.PCA(f)
  }

  caPCA
}
#' List files with filtered reaction coordinates
.list.reactionCoords <- function() {
  list.files(pattern = "reaction_coords_[[:digit:]]+.coords$")
}

.check.projectPath <- function() {
  if (exists("path", envir=.project_cache, inherits=FALSE)) {
    path <- get("path", envir=.project_cache)
    if (getwd() != path) {
      warning(.warnings$diverging_path)
      setwd(path)
    }
  } else {
    stop(.warnings$no_init)
  }
  NULL
}

#' create a README file
#'
#' Creates a README file to describe all automatically
#' generated files in the project path.
.make.readme <- function() {
  sec <- function(title, level=1) {
    separator <- paste(rep("\n", 1), collapse="")
    indicator <- paste(rep("#", level), collapse="")
    paste(separator, paste(indicator, title), sep="")
  }
  txt <- function(text) {
    paste(text)
  }

  pd <- project()

  proj_summary <- ""

  add_doc <- function(elem, doc) {
    if ((! is.null(elem)) & (! length(elem) == 0)) {
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
  if ( ! is.null(pd$caDists) & length(pd$caDists) > 0) {
    proj_summary <- c(proj_summary, sec("C$_\\alpha$-distances"))
    for (d in pd$caDists) {
      add_doc(d, txt(d))
    }
  }
  if ( ! is.null(pd$caPCA) & length(pd$caPCA) > 0) {
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
      proj_summary <<- c(proj_summary,
                         sec(fname, level="2"),
                         txt(readLines(paste(strsplit(fname, "\\.")[[1]][1],
                                             "desc",
                                             sep="."))))
    })
  }


  cat(paste(c(proj_summary, ""), sep="", collapse="\n"),
      file="README.Rmd")
  opts <- list()
  opts["toc"] <- TRUE
  sink("/dev/null")
  suppressMessages(rmarkdown::render("README.Rmd",
                                     "pdf_document",
                                     output_options = opts))
  sink()
}

#' update project information
.update <- function() {
  .check.projectPath()

  project <- list()
  project$ref <- .list.reference()
  project$traj <- .list.traj()
  project$dihedrals <- .list.dihedrals()
  project$dPCAplus <- .list.dPCAplus()
  project$caDists <- .list.caDists()
  project$caPCA <- .list.caPCA()
  project$reactionCoords <- .list.reactionCoords()

  assign("project", project, envir=.project_cache)

  .make.readme()
}

#' initialize project
#'
#' Switches to given path and initializes prodyna-project.
#' Existing files are detected and any new files generated by
#' prodyna will be created under this directory.
#'
#' @export
init <- function(path=".") {
  assign("path",
         normalizePath(path),
         envir=.project_cache)
  .update()
}
