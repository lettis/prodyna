
#
# TODO: package documentation & license(-file)
#

#' warning messages
.warnings <- list("no_init" =
                    "No project has been initialized so far. Call prodyna::initialize(...)",
                  "diverging_path" =
                    "Current directory and project paths diverge. Resetting current directory to project path")


#### handling external binaries

.default_binaries <- list("gmx"="/usr/local/gromacs/bin/gmx",
                          "fastpca"="/usr/local/bin/fastpca",
                          "clustering"="/usr/local/bin/clustering",
                          "bash"="/bin/bash",
                          "awk"="/usr/bin/awk",
                          "paste"="/usr/bin/paste")

.binaries <- new.env()

#' Set path to binary
#'
#' @param name Name of the binary
#' @param path Absolute system path to binary
#' @examples
#' set.binary("bash", "/bin/bash")
#' @export
set.binary <- function(name, path) {
  assign(name, path, envir=.binaries)
}

#' Get path to binary
#'
#' @param name Name of the binary
#' @export
get.binary <- function(name) {
  if (exists(name, envir=.binaries, inherits=FALSE)) {
    get(name, envir=.binaries)
  } else {
    .default_binaries[[name]]
  }
}

#' List names and paths to binaries called by prodyna
#' @export
list.binaries <- function() {
  lapply(names(.default_binaries), function(n) {
    c(n, get.binary(n))
  })
}


#### project management

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
    if ( ! is.null(elem)) {
      proj_summary <<- c(proj_summary, doc)
    }
  }

  add_doc(pd$ref, c(sec("Reference Structure"),
                    txt(pd$ref)))
  add_doc(pd$traj, c(sec("Trajectory"),
                     txt(pd$traj)))
  add_doc(pd$dihedrals, c(sec("Dihedrals"),
                          txt(pd$dihedrals)))
  add_doc(pd$dPCAplus, sec("dPCA+ Analysis"))
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
  if ( ! is.null(pd$caDists)) {
    proj_summary <- c(proj_summary, sec("C$_\\alpha$-distances"))
    for (d in pd$caDists) {
      add_doc(d, txt(d))
    }
  }
  if ( ! is.null(pd$caPCA)) {
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
initialize <- function(path=".") {
  assign("path",
         normalizePath(path),
         envir=.project_cache)
  .update()
}



#### generator functions

#' generate dihedrals, return updated project information
#' @export
generate.dihedrals <- function() {
  .check.projectPath()
  # get project description
  pd <- project()
  filter.backbone <- function(pdb) {
    pdb$atom$elety == "N" |
      pdb$atom$elety == "CA" |
      pdb$atom$elety == "C"
  }
  pdb_ref <- bio3d::read.pdb(pd$ref[[1]])
  # correct for custom residue types
  pdb_ref$calpha <- pdb_ref$atom$elety == "CA"
  n_atoms <- length(pdb_ref$calpha)
  calpha_indices <- pdb_ref$atom$eleno[pdb_ref$calpha]
  # remove first and last calpha
  # (incomplete dihedrals at terminal regions)
  calpha_indices <- head(calpha_indices, -1)[-1]
  bb <- pdb_ref$atom$eleno[filter.backbone(pdb_ref)]
  dih_indices <- lapply(calpha_indices, function(ca){
    i <- which(bb == ca)
    phi <- c(bb[i-2], bb[i-1], bb[i], bb[i+1])
    psi <- c(bb[i-1], bb[i], bb[i+1], bb[i+2])
    c(phi, psi)
  })
  dih_indices <- do.call("c", dih_indices)
  # prepare index file to identify backbone-dihedral atoms
  tmp_ndx <- "tmp_phipsi.ndx"
  unlink(tmp_ndx)
  write("[PhiPsi]", tmp_ndx)
  write(dih_indices, tmp_ndx, append=TRUE, ncolumns=4)
  # compute dihedrals
  fname_dihedrals <- paste(pd$traj, ".dih", sep="")
  fname_xvg <- paste(fname_dihedrals, ".xvg", sep="")
  unlink(fname_xvg)
  print("running GMX to generate dihedrals")
  system2(gmx.binary, c("angle",
                        " -f ",
                        pd$traj,
                        " -n ",
                        tmp_ndx,
                        " -ov ",
                        fname_xvg,
                        " -type dihedral -all"),
          stdout = TRUE,
          stderr = TRUE)
  # cleanup
  dih <- fread(fname_xvg, drop=c(1,2))
  fwrite(dih,
         fname_dihedrals,
         sep=" ",
         quote=FALSE,
         row.names=FALSE,
         col.names=FALSE)
  unlink(tmp_ndx)
  unlink("angdist.xvg")
  unlink(fname_xvg)
  print("done")

  pd$dihedrals <- list.dihedrals()

  .update()
}

#' generate trajectory with C\eqn{\alpha} distances as coordinates
#' @export
generate.caDistances <- function(residue.mindist=4, residue.maxdist=NULL) {
  .check.projectPath()
  # get project description
  pd <- project()

  pdb <- bio3d::read.pdb(pd$ref)
  calpha_ndx <- pdb$atom$eleno[pdb$atom$elety == "CA"]
  n_res <- length(calpha_ndx)
  res_ndx <- 1:n_res
  if (is.null(residue.maxdist)) {
    residue.maxdist <- n_res
  }
  ca_comb <- combn(calpha_ndx, 2)
  ca_pairs <- do.call("rbind",
                      Filter(Negate(is.null),
                             lapply(1:ncol(ca_comb), function(i) {
                               res1 <- res_ndx[calpha_ndx == ca_comb[1,i]]
                               res2 <- res_ndx[calpha_ndx == ca_comb[2,i]]
                               pair <- NULL
                               if (res1 < res2) {
                                 d_res <- res2-res1
                                 if (residue.mindist <= d_res & d_res <= residue.maxdist) {
                                   pair <- c(ca_comb[1,i], ca_comb[2,i])
                                 }
                               }
                               pair
                             })))

  ca_dist_fname <- caDistFilename(pd$traj,
                                  residue.mindist,
                                  residue.maxdist)
  ndx_fname <- paste(ca_dist_fname,
                     ".ndx",
                     sep="")
  # write new index file with distance pairs
  unlink(ndx_fname)
  for (i in 1:nrow(ca_pairs)) {
    write(paste("[ Dist",
                i,
                " ]\n",
                ca_pairs[i,1],
                " ",
                ca_pairs[i,2],
                sep=""),
          file=ndx_fname,
          append=TRUE)
  }
  # compute distances with gromacs
  system2(gmx.binary, c("distance",
                        "-f",
                        pd$traj,
                        "-n",
                        ndx_fname,
                        "-oall",
                        "-select",
                        seq(0, (nrow(ca_pairs)-1))))
  # reformat data
  system(paste("grep -v \"#\" dist.xvg | grep -v \"@\" ",
               "| ",
               awk.binary,
               " '{for(i=2; i <=NF; ++i)",
               "{printf(\" %s\", $i)} printf(\"\\n\")}'",
               " > ",
               ca_dist_fname,
               sep=""))
  # remove intermediate file
  unlink("dist.xvg")
  # update project information
  if (is.null(pd$caDists)) {
    pd$caDists <- ca_dist_fname
  } else {
    pd$caDists <- c(pd$caDists, ca_dist_fname)
  }

  .update()
}

#' run dPCA+
#' @export
run.dPCAplus <- function(corr=FALSE) {
  .check.projectPath()
  # get project information
  pd <- project()
  if (!("dihedrals" %in% names(pd)) | length(pd$dihedrals) == 0) {
    pd <- generate.dihedrals(pd)
  }
  pd$dPCAplus <- list()
  if (corr) {
    postfix <- "n"
  } else {
    postfix <- ""
  }
  for (name in c("proj", "vec", "val", "cov", "stats")) {
    pd$dPCAplus[[paste(name, postfix, sep="")]] <- paste(pd$dihedrals,
                                                         ".",
                                                         name,
                                                         postfix,
                                                         sep="")
  }
  # run PCA
  if (corr) {
    params <- c("-f",
                pd$dihedrals,
                "-p",
                pd$dPCAplus$projn,
                "-c",
                pd$dPCAplus$covn,
                "-v",
                pd$dPCAplus$vecn,
                "-l",
                pd$dPCAplus$valn,
                "-s",
                pd$dPCAplus$statsn,
                "-P",
                "-N")
  } else {
    params <- c("-f",
                pd$dihedrals,
                "-p",
                pd$dPCAplus$proj,
                "-c",
                pd$dPCAplus$cov,
                "-v",
                pd$dPCAplus$vec,
                "-l",
                pd$dPCAplus$val,
                "-s",
                pd$dPCAplus$stats,
                "-P")
  }
  # run PCA
  system2(paste(fastpca.dir, "fastpca", sep=""), args=params)

  .update()
}

#' run C\eqn{\alpha} distance PCA
#' @export
run.caPCA <- function(residue.mindist=4, residue.maxdist=NULL, corr=FALSE) {
  .check.projectPath()
  # get project information
  pd <- project()
  ca_dist_fname <- caDistFilename(pd$traj,
                                  residue.mindist,
                                  min(residue.maxdist,
                                      get.nResidues(pd$ref)))
  if (! ca_dist_fname %in% pd$caDists) {
    pd <- generate.caDistances(pd,
                               residue.mindist=residue.mindist,
                               residue.maxdist=residue.maxdist)
  }
  if (is.null(pd$caPCA)) {
    pd$caPCA <- list()
  }
  pca_label <- paste(residue.mindist,
                     "_",
                     residue.maxdist,
                     sep="")
  pd$caPCA[[pca_label]] <- list()
  if (corr) {
    postfix <- "n"
  } else {
    postfix <- ""
  }
  for (name in c("proj", "vec", "val", "cov", "stats")) {
    pd$caPCA[[pca_label]][[paste(name, postfix, sep="")]] <- paste(ca_dist_fname,
                                                                   ".",
                                                                   name,
                                                                   postfix,
                                                                   sep="")
  }
  # run PCA
  if (corr) {
    params <- c("-f",
                ca_dist_fname,
                "-p",
                pd$caPCA[[pca_label]]$projn,
                "-c",
                pd$caPCA[[pca_label]]$covn,
                "-v",
                pd$caPCA[[pca_label]]$vecn,
                "-l",
                pd$caPCA[[pca_label]]$valn,
                "-s",
                pd$caPCA[[pca_label]]$statsn,
                "-N")
  } else {
    params <- c("-f",
                ca_dist_fname,
                "-p",
                pd$caPCA[[pca_label]]$proj,
                "-c",
                pd$caPCA[[pca_label]]$cov,
                "-v",
                pd$caPCA[[pca_label]]$vec,
                "-l",
                pd$caPCA[[pca_label]]$val,
                "-s",
                pd$caPCA[[pca_label]]$stats)
  }
  # run PCA
  system2(fastpca.binary, args=params)

  .update()
}

#' select subspace projection
#' @export
select.reactionCoordinates <- function(coords, columns) {
  .check.projectPath()
  if (xor(is.list(coords), is.list(columns))) {
    stop("either coords and columns are both lists (of same length), or none is")
  }
  if ( ! is.list(coords)) {
    coords <- list(coords)
    columns <- list(columns)
  }
  # get project information
  pd <- project()
  # generate filenames
  id <- length(pd$reactionCoords) + 1
  bname <- "reaction_coords_"
  desc_fname <- paste(bname, id, ".desc", sep="")
  coords_fname <- paste(bname, id, ".coords", sep="")

  # generate description
  desc <- do.call(c, lapply(1:length(coords),
                            function(i) {
                              paste(coords[[i]],
                                    paste(columns[[i]],
                                          collapse=" "))
                            }))
  cat(desc, file=desc_fname, sep="\n")

  # use awk, paste and bash to filter reaction coords
  pipes <- do.call(c, lapply(1:length(coords), function(i) {
    cols <- paste(paste("\\$", columns[[i]], sep=""), collapse=", ")
    paste("<(",
          awk.binary,
          "'{print",
          cols,
          "}'",
          coords[[i]],
          ")")
  }))
  cmd <- paste(paste(bash.binary, " -c \"", sep=""),
               paste.binary,
               "-d ' '",
               paste(pipes, collapse=" "),
               ">",
               coords_fname,
               "\"")
  system(cmd)

  .update()
}

#### file IO

#' read dihedrals
#' @return data.frame with \eqn{\phi} and \eqn{\psi} angles
#' @param resno Residue number. If NULL: read dihedrals from all residues.
#' @export
read.dihedrals <- function(resno=NULL) {
  .check.projectPath()
  # get project information
  pd <- project()
  if (is.null(resno)) {
    dih <- data.table::fread(pd$dihedrals)
    n_dih <- dim(dih)[2]/2
    res_ndx <- do.call("c", lapply(2:(n_dih+1), function(i){rep(i,2)}))
  } else {
    # TODO resno as vector?
    dih <- data.table::fread(pd$dihedrals, select=c(2*(resno-1)-1, 2*(resno-1)))
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
#' p <- prodyna::project()
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

#### plotting

#' plot Ramachandran plot
#' @param resno Residue number.
#' @export
plot.ramachandran <- function(resno) {
  .check.projectPath()
  pd <- project()
  suppressMessages(require(ggplot2))
  dih <- read.dihedrals(pd, resno)
  phi <- dih[[paste("phi", resno, sep="")]]
  psi <- dih[[paste("psi", resno, sep="")]]
  p <- ggplot(data.frame(phi, psi)) +
    stat_bin2d(aes(x=phi, y=psi), bins=180) +
    scale_fill_distiller(palette="YlGnBu",
                         trans=reverselog_trans()) +
    xlim(-180,180) +
    ylim(-180,180) +
    theme(legend.position = "none")

  p
}

#' plot 2d-proj, 1d-proj and eigenvector content for given PCA
#' @export
plot.pcaOverview <- function(pca, pcs, corr=FALSE) {
  .check.projectPath()
  suppressMessages(require(data.table))
  suppressMessages(require(GGally))
  suppressMessages(require(ggplot2))
  if (corr) {
    check.filesExist(c(pca$projn, pca$vecn))
    proj <- fread(pca$projn, select=pcs)
    vecs <- fread(pca$vecn, select=pcs)
  } else {
    check.filesExist(c(pca$proj, pca$vec))
    proj <- fread(pca$proj, select=pcs)
    vecs <- fread(pca$vec, select=pcs)
  }
  vec_names <- names(vecs)
  logproj <-
    n_dih <- dim(vecs)[1]
  vecs <- do.call("rbind",
                  lapply(vec_names, function(n) {
                    data.frame(dih=vecs[[n]], v=n, ndx=1:n_dih)
                  }))
  plt <- ggpairs(proj,
                 lower="blank",
                 upper="blank")
  seq <- 1:ncol(proj)
  for (x in seq){
    for (y in seq) {
      if (x == y) {
        # diag
        plt <- putPlot(plt,
                       ggplot(proj,
                              aes_string(x=names(proj)[x])) +
                         stat_bin(bins=200,
                                  aes(y=-log(..density..)),
                                  geom="line",
                                  position="identity"),
                       x,
                       x)
      }
      if (x < y) {
        # lower
        plt <- putPlot(plt,
                       ggplot(proj,
                              aes_string(x=names(proj)[x],
                                         y=names(proj)[y])) +
                         stat_bin2d(bins=200) +
                         scale_fill_distiller(palette="YlGnBu",
                                              trans=reverselog_trans()),
                       y,
                       x)
      } else if (y < x) {
        # upper
        plt <- putPlot(plt,
                       ggplot(vecs[vecs$v==vec_names[x] | vecs$v==vec_names[y],],
                              aes(x=ndx,
                                  y=abs(dih),
                                  group=v,
                                  color=v)) +
                         geom_line(size=2),
                       y,
                       x)
      }
    }
  }

  plt
}

#' plot cumulative fluctuations
#' Takes all available PCAs and plots their cumulative fluctuations in a single plot.
#' @export
plot.cumFlucts <- function() {
  .check.projectPath()
  pd <- project()

  #TODO different caPCAs!!!

  # check different kinds of PCA for results
  cfs_labels <- c("dPCA+",
                  "dPCA+ (corr)",
                  "caPCA",
                  "caPCA (corr)")
  cfs_filenames <- c(pd$dPCAplus$val,
                     pd$dPCAplus$valn,
                     pd$caPCA$val,
                     pd$caPCA$valn)
  max_length <- 0
  cfs <- list()
  cfs_labels_selected <- list()
  for (i in 1:length(cfs_labels)) {
    f <- cfs_filenames[i]
    if (file.exists(f)) {
      cf <- fread(f)$V1
      cf <- cumsum(cf/sum(cf))
      max_length <- max(max_length, length(cf))
      cfs[[i]] <- data.frame(1:max_length, cf, cfs_labels[[i]])
      colnames(cfs[[i]]) <- c("PC", "cumfluct", "method")
    }
  }
  for (i in length(cfs)) {
    cfs[[i]] <- cfs[[i]][1:max_length,]
  }
  df <- do.call("rbind", cfs)
  # plot result
  ggplot(df) +
    geom_line(aes(x=PC, y=cumfluct, color=method),
              size=2) +
    xlim(1, max_length)
}

## plot autocorrelation of observables
##   coords: input file with coordinates
##   columns: vector of column indices to use
##   corrlength: extend of autocorrelation computation.
##               if < 1, ratio of number of data points
plot.autocorr <- function(coords, columns, corrlength=0.25) {
  #TODO: finish
}
