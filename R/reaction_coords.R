
#' generate dihedrals, return updated project information
#' @param skipCA Vector of CA indices to be skipped. Default: first and last.
#' @param ignoreCache Ignore cached files and recalculate in any case.
#' @export
generate.dihedrals <- function(skipCA=NULL, ignoreCache=FALSE) {
  .check.projectPath()
  # get project description
  pd <- project()

  if (length(pd$dihedrals) == 0 | ignoreCache) {
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
    if (is.null(skipCA)) {
      calpha_indices <- head(calpha_indices, -1)[-1]
    } else {
      calpha_indices <- calpha_indices[!(1:length(calpha_indices) %in% skipCA)]
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
    fname_dihedrals <- paste(pd$traj, ".dih", sep="")
    fname_xvg <- paste(fname_dihedrals, ".xvg", sep="")
    unlink(fname_xvg)
    print("running GMX to generate dihedrals")
    system2(get.binary("gmx"), c("angle",
                                 " -f ",
                                 pd$traj,
                                 " -n ",
                                 tmp_ndx,
                                 " -ov ",
                                 fname_xvg,
                                 " -type dihedral -all"),
            stdout = TRUE,
            stderr = TRUE)
    # streamed xvg -> dih conversion
    cmd <- paste(get.binary("awk"),
                 "'!/^#.*/ && !/^@.*/ {for(i=3; i<=NF; ++i)",
                 #"'{for(i=3; i<=NF; ++i)",
                 "printf(\" %s\", $i); printf(\"\\n\")}'",
                 fname_xvg,
                 ">",
                 fname_dihedrals)
    system(cmd)
    # cleanup
    unlink(tmp_ndx)
    unlink("angdist.xvg")
    unlink(fname_xvg)
    print("done")
  }

  .update()
}

#' run dPCA+
#'
#' Runs a dPCA+ analysis for the given project.
#' If dihedrals have not been generated yet, they will be automatically
#' using the default settings (remove first and last C\eqn{\alpha} to get
#' rid of the end-caps).
#' @param corr Use correlation instead of covariance. Effectively whitens the data before doing analysis.
#' @param ignoreCache Ignore cached files and recalculate in any case.
#' @export
run.dPCAplus <- function(corr=FALSE, ignoreCache=FALSE) {
  .check.projectPath()
  # get project information
  pd <- project()
  if (!("dihedrals" %in% names(pd)) | length(pd$dihedrals) == 0) {
    generate.dihedrals()
  }

  dpca_keys <- names(pd$dPCAplus)
  missing_files <- any(do.call(c, lapply(dpca_keys[grep("*n$",
                                                        dpca_keys,
                                                        invert=(!corr))],
                                         function(x) {length(pd$dPCAplus[[x]]) == 0})))
  if (is.null(pd$dPCAplus) |
      missing_files |
      ignoreCache) {
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
    system2(get.binary("fastpca"), args=params)
  }

  .update()
}



#######
# TODO: caching for all functions below!
#######

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
