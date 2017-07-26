
#' Compute dihedral angles.
#'
#' Dihedral angles for the given trajectory are computed using gmx.
#' The result is saved to <traj>.dih. A file <traj>.dih.info is generated
#' containing additional information.
#'
#' @param ref Character, name of the PDB file describing the reference structure.
#' @param traj Character, name of the XTC file describing the trajectory.
#' @param skipCA Numerical vector, C\eqn{_\alpha} indices to be skipped.
#'   If \code{NULL} (default) first and last C\eqn{_\alpha} indices are skipped
#'   (incomplete dihedrals at terminal regions).
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output file
#'   already exists (default: \code{FALSE}). Generates a warning if cached files
#'   exist.
#' @export
generate.dihedrals <- function(ref, traj, skipCA=NULL, ignoreCache=FALSE) {
  # output file
  fname_dihedrals <- paste(traj, ".dih", sep="")
  fname_dihedrals_info <- paste(traj, ".dih.info", sep="")

  #if (file.exists(fname_dihedrals) & file.exists(fname_dihedrals_info)) {
  #  dih_info <- read.dihedrals.info(fname_dihedrals_info)

  if (!file.exists(fname_dihedrals)|ignoreCache) {
    filter.backbone <- function(pdb) {
      pdb$atom$elety == "N" |
      pdb$atom$elety == "CA" |
      pdb$atom$elety == "C"
    }
    pdb_ref <- bio3d::read.pdb(ref)

    # correct for custom residue types
    pdb_ref$calpha <- pdb_ref$atom$elety == "CA"
    n_atoms <- length(pdb_ref$calpha)
    calpha_indices <- pdb_ref$atom$eleno[pdb_ref$calpha]

    # remove skipCA indices
    if (is.null(skipCA)) {
      # default: remove first and last C_alpha
      resnos <- 2:(length(calpha_indices)-1)
      calpha_indices <- head(calpha_indices, -1)[-1]
    } else {
      resnos <- (1:length(calpha_indices))[!(1:length(calpha_indices) %in% skipCA)]
      calpha_indices <- calpha_indices[resnos]
    }
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
    fname_xvg <- paste(fname_dihedrals, ".xvg", sep="")
    unlink(fname_xvg)

    print("Running GMX to generate dihedrals.")

    run.cmd("gmx", args = c("angle",
                            " -f ",
                            traj,
                            " -n ",
                            tmp_ndx,
                            " -ov ",
                            fname_xvg,
                            " -type dihedral -all"))

    # streamed xvg -> dih conversion
    cmds <- paste(get.binary("awk"),
                 "'!/^#.*/ && !/^@.*/ {for(i=3; i<=NF; ++i)",
                 #"'{for(i=3; i<=NF; ++i)",
                 "printf(\" %s\", $i); printf(\"\\n\")}'",
                 fname_xvg,
                 ">",
                 fname_dihedrals)
    run.cmds(cmds)

    # write .dih.info file
    write(paste(paste("reference",  ref),
                paste("trajectory", traj),
                paste("residuals", paste(resnos, collapse=" ")),
                paste("skippedDefault", is.null(skipCA)),
                sep = "\n"
                ),
          file = fname_dihedrals_info
          )

    # cleanup
    unlink(tmp_ndx)
    unlink("angdist.xvg")
    unlink(fname_xvg)
    print("Done.")
  } else {
    warning(msg.warn("caching", arg="generate.dihedrals"))
  }
}

#' Generate cos/sin transformed dihedrals.
#'
#' @param dihedrals Character, name of .dih file containing dihedral angles.
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output file
#'   already exists (default: \code{FALSE}).
#' @export
generate.cossinTransform <- function(dihedrals, ignoreCache=FALSE) {
  output <- paste(dihedrals, "cossin", sep=".")
  if (!file.exists(output)|ignoreCache) {
    streamed_cossin_transform(dihedrals, output)
  }
}


#' PCA
#'
#' Run PCA on given coordinates.
#'
#' @param coords Character, name of file containing the coordinates.
#' @param corr Logical, if \code{TRUE} use correlation instead of covariance.
#'   Effectively whitens the data before doing analysis (default:
#'   \code{FALSE}).
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output files
#'   already exist (default: \code{FALSE}).
#' @param additionalParams Character vector, additional parameters for 'fastpca'
#'   (default: \code{NULL}).
#' @export
run.PCA <- function(coords, corr=FALSE, ignoreCache=FALSE, additionalParams=NULL) {
  results <- c("proj", "vec", "val", "cov", "stats")
  if (corr) {
    results <- paste(results, "n", sep="")
  }
  results <- paste(coords, results, sep=".")
  # results already computed?
  cached_results_missing <- (!all(file.exists(results)))

  if (ignoreCache | cached_results_missing) {
    # setup 'fastpca' parameters
    params <- c("-f", "-p", "-v", "-l", "-c", "-s")
    params <- c(rbind(params, c(coords, results)))
    if (corr) {
      params <- c(params, "-N")
    }
    params <- c(params, additionalParams)
    print("Running 'fastpca'..")
    run.cmd("fastpca", args=params)
    print(".. done.")
  } else {
    print("PCA output files already exist.")
  }
}


#' dPCA+
#'
#' Runs a dPCA+ analysis on dihedral angles.
#'
#' @param dihedrals Character, name of .dih file containing dihedral angles.
#' @param corr Logical, if \code{TRUE} use correlation instead of covariance.
#'  Effectively whitens the data before doing analysis (default: \code{FALSE}).
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output files
#'   already exists (default: \code{FALSE}).
#' @export
run.dPCAplus <- function(dihedrals, corr=FALSE, ignoreCache=FALSE) {
  run.PCA(coords=dihedrals,
          corr=corr,
          ignoreCache=ignoreCache,
          additionalParams = c("-P"))
}


### TODO: caching for all functions below!


#' generate trajectory with C\eqn{\alpha} distances as coordinates
#' @export
generate.caDistances <- function(residue.mindist=4, residue.maxdist=NULL) {
  .check.projectPath()
  # get project description
  pd <- projectInfo()

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

#' run C\eqn{\alpha} distance PCA
#' @export
run.caPCA <- function(residue.mindist=4, residue.maxdist=NULL, corr=FALSE) {
  .check.projectPath()
  # get project information
  pd <- projectInfo()
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

#' Generate reaction coordinates.
#'
#' Select subspace projection.
#'
#' @param coords Either filename or list of filenames of coordinates.
#' @param columns Select columns. Either single vector or
#'                list of vectors as selection per coord file.
#' @param output Output filename (default: NULL). If NULL, filename will be generated.
#' @export
generate.reactionCoordinates <- function(coords, columns, output=NULL) {
  if (xor(is.list(coords), is.list(columns))) {
    stop("either coords and columns are both lists (of same length), or none is")
  }
  if ( ! is.list(coords)) {
    coords <- list(coords)
    columns <- list(columns)
  }
  if (is.null(output)) {
    # generate filenames
    id <- length(pd$reactionCoords) + 1
    output <- get.fullPath(paste("reaction_coords_", id, sep=""))
    coords_fname <- paste(output, ".coords", sep="")
  } else {
    coords_fname <- output
  }

  desc_fname <- paste(output, ".desc", sep="")

  if (file.exists(desc_fname)|file.exists(coords_fname)) {
    warning("no reaction coordinates generated: file exists")
  } else {
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
            get.binary("awk"),
            "'{print",
            cols,
            "}'",
            coords[[i]],
            ")")
    }))
    cmd <- paste(paste(get.binary("bash"), " -c \"", sep=""),
                 get.binary("paste"),
                 "-d ' '",
                 paste(pipes, collapse=" "),
                 ">",
                 coords_fname,
                 "\"")
    system(cmd)
  }
}

#' filter a data set
#'
#' @param coords Filename of coordinates (ASCII data)
#' @param filter Some integer N. If 'states' is NULL (default),
#'               every N-th frame will be selected.
#'               If 'states' points to file, N is selected state.
#' @param states Name of state trajectory. default: NULL, i.e. not used.
#' @param output Filename of output. Default: NULL,
#'               i.e. construct a name from input parameters.
#' @param ignoreCache Recompute even if results already exist.
#' @export
generate.filteredCoordinates <- function(coords, filter, states=NULL, output=NULL, ignoreCache=FALSE) {
  if ( ! file.exists(coords)) {
    coords <- get.fullPath(coords)
  }
  if (is.null(states)) {
    if (is.null(output)) {
      output <- paste(coords, ".every", filter, sep="")
    }
    cmd <- paste(get.binary("awk"),
                 " 'NR % ",
                 filter,
                 " == 0 {print $0}' ",
                 coords,
                 " > ",
                 output,
                 sep="")
  } else {
    if (is.null(output)) {
      output <- paste(coords, ".state_", filter, sep="")
    }
    cmd <- paste(get.binary("clustering"),
                 " filter -s ",
                 states,
                 " -c ",
                 coords,
                 " -o ",
                 output,
                 " -S ",
                 filter,
                 sep="")
  }
  if ( ( ! file.exists(output)) | ignoreCache) {
    # run filter function
    system(cmd)
  }
}
