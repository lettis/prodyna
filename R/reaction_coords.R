
#' Compute dihedral angles.
#'
#' Dihedral angles for the given trajectory are computed using gmx.
#' The result is saved to <traj>.dih. A file <traj>.dih.info is generated
#' containing additional information (see \code{\link{write.dihedrals.info}}).
#'
#' If \code{ignoreCache} is set to \code{FALSE}, both output files already exist
#' and the skipped residuals in the .dih.info match \code{skipCA}, no
#' computation is performed.
#' In any other case, the results are recomputed (and potentially existing
#' output files are overwritten).
#'
#' @param ref Character, name of the PDB file describing the reference structure.
#' @param traj Character, name of the XTC file describing the trajectory.
#' @param skipCA Numeric vector, C\eqn{_\alpha} indices to be skipped in
#'   increasing order. \cr
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

  # check if cached files exist
  if (!ignoreCache &&
      file.exists(fname_dihedrals) &&
      file.exists(fname_dihedrals_info)) {

    dihInfo <- read.dihedrals.info(fname_dihedrals_info)
    nRes    <- dihInfo$nRes
    skipped <- dihInfo$skippedResnos

    # check if the same residuals should be skipped
    if ((is.null(skipCA) &&
         identical(skipped, c(1, nRes))) ||
        (!is.null(skipCA) &&
         identical(skipped, skipCA) )) {

      return(warning(msg("caching", arg="generate.dihedrals")))
    }
  }

  # if reference is .gro file, convert to .pdb first
  if (substring(ref, nchar(ref)-2, nchar(ref)) == "gro") {

    ref_pdb <- paste(substring(ref, 1, nchar(ref)-4), "pdb", sep=".")

    run.cmd(cmd = get.binary("gmx"),
            args = c("editconf",
                     "-f",
                     ref,
                     "-o",
                     ref_pdb))
    ref <- ref_pdb
  }

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

  nRes <- length(calpha_indices)


  if (is.null(skipCA)) {
    skipCA <- c(1, nRes)               # default: remove first and last C_alpha
  } else {
    skipCA <- skipCA[skipCA <= nRes]   # remove out of bound indices
  }

  # remove skipCA indices
  resnos <- (1:nRes)[-skipCA]
  calpha_indices <- calpha_indices[resnos]

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

  # in case of an error delete all generated files
  onError <- function() {unlink(c(fname_dihedrals,
                                  fname_dihedrals_info,
                                  fname_xvg,
                                  tmp_ndx,
                                  "angdist.xvg"))}

  message("Running GMX to generate dihedrals.")

  run.cmd(get.binary("gmx"),
          args = c("angle",
                   " -f ",
                   traj,
                   " -n ",
                   tmp_ndx,
                   " -ov ",
                   fname_xvg,
                   " -type dihedral -all"),
          onError= onError)

  # streamed xvg -> dih conversion
  cmds <- paste(get.binary("awk"),
               "'!/^#.*/ && !/^@.*/ {for(i=3; i<=NF; ++i)",
               #"'{for(i=3; i<=NF; ++i)",
               "printf(\" %s\", $i); printf(\"\\n\")}'",
               fname_xvg,
               ">",
               fname_dihedrals)
  run.cmds(cmds, onError=onError)

  # write .dih.info file
  write.dihedrals.info(ref, traj, nRes, skipCA, fname_dihedrals_info)

  # cleanup
  unlink(tmp_ndx)
  unlink("angdist.xvg")
  unlink(fname_xvg)
  message("Done.")
}

#' Cosine/Sine Transformation.
#'
#' Generate cos/sin-transformed dihedrals.
#'
#' @param dihedrals Character, name of .dih file containing dihedral angles.
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output file
#'   already exists (default: \code{FALSE}).
#' @export
generate.cossinTransform <- function(dihedrals, ignoreCache=FALSE) {

  output <- paste(dihedrals, "cossin", sep=".")

  if (!ignoreCache && file.exists(output)) {
    warning(msg("caching", arg="generate.cossinTransform"))
  } else {
    if(!file.exists(dihedrals)) {
      stop(msg("missingFile", dihedrals))
    } else {
      streamed_cossin_transform(dihedrals, output)
    }
  }
}


#' PCA
#'
#' Run PCA on the given coordinates using fastpca. \cr
#' The following files are created:
#' \itemize{
#' \item if \code{corr} is \code{FALSE}: <coords>.proj, <coords>.vec,
#'       <coords>.val, <coords>.cov, <coords>.stats
#' \item if \code{corr} is \code{TRUE}: <coords>.projn, <coords>.vecn,
#'       <coords>.valn, <coords>.covn, <coords>.statsn
#' }
#'
#' @param coords Character, name of file containing the coordinates.
#' @param corr Logical, if \code{TRUE} use correlation instead of covariance.
#'   Effectively whitens the data before doing the analysis (default:
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
  outputFilesExist <- all(file.exists(results))

  if (!ignoreCache && outputFilesExist) {
    warning(msg("caching", "run.PCA"))
  } else {
    if(!file.exists(coords)) {
      stop(msg("missingFile", coords))
    }
    # setup 'fastpca' parameters
    params <- c("-f", "-p", "-v", "-l", "-c", "-s")
    params <- c(rbind(params, c(coords, results)))
    if (corr) {
      params <- c(params, "-N")
    }
    params <- c(params, additionalParams)
    message("Running 'fastpca'... ", appendLF=F)

    run.cmd(get.binary("fastpca"),
            args    = params,
            onError = function(){unlink(results)})

    message("done.")
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
#'   already exist (default: \code{FALSE}).
#' @export
run.dPCAplus <- function(dihedrals, corr=FALSE, ignoreCache=FALSE) {
  run.PCA(coords=dihedrals,
          corr=corr,
          ignoreCache=ignoreCache,
          additionalParams = c("-P"))
}


#' C\eqn{\alpha} distances.
#'
#' Generate a trajectory with C\eqn{\alpha} distances as coordinates.
#'
#' @param ref Character, name of the PDB file describing the reference structure.
#' @param traj Character, name of the XTC file describing the trajectory.
#' @param residue.mindist (default: 4)
#' @param residue.maxdist (default:\code{NULL})
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output file
#'   already exists (default: \code{FALSE}).
#' @importFrom bio3d read.pdb
#' @export
generate.caDistances <- function(ref, traj, residue.mindist=4, residue.maxdist=NULL, ignoreCache=FALSE) {

  pdb <- read.pdb(ref)
  calpha_ndx <- pdb$atom$eleno[pdb$atom$elety == "CA"]
  n_res <- length(calpha_ndx)

  if (is.null(residue.maxdist)) {
    residue.maxdist <- n_res
  }

  fname_ca_dists <- paste(traj, "caDist", sep=".")
  fname_ca_dists_info <- paste(fname_ca_dists, "info", sep=".")

  # check if cached files exist
  if (!ignoreCache &&
      file.exists(fname_ca_dists) &&
      file.exists(fname_ca_dists_info)) {

    caDistInfo <- read.caDistances.info(fname_ca_dists_info)

    # check if parameters are the same
    if (caDistInfo$mindist == residue.mindist &&
        caDistInfo$maxdist == residue.maxdist) {

      return(warning(msg("caching", arg="generate.caDistances")))
    }
  }

  # filter all index pairs (i,j), i, j in {1,.., n_res} with i < j
  filter <- function(d_res) {residue.mindist <= d_res && d_res <= residue.maxdist}

  ca_pairs <- Reduce(
               function(l1, i) {
                 l2 <-Reduce(function(l,j) {
                               d_res <- j-i
                               if (filter(d_res)) {
                                 return(c(l, list(c(calpha_ndx[i], calpha_ndx[j]))))
                               } else {
                                 return(l)}
                                },
                             x    = (i+1):n_res,
                             init = list())
                 return(c(l1, l2))},
               x    = 1:(n_res-1),
               init = list())

  ca_pairs <- do.call("rbind", ca_pairs)

  fname_ndx <- paste(fname_ca_dists, ".ndx", sep="")

  # write new index file with distance pairs
  unlink(fname_ndx)
  for (i in 1:nrow(ca_pairs)) {
    write(paste("[ Dist",
                i,
                " ]\n",
                ca_pairs[i,1],
                " ",
                ca_pairs[i,2],
                sep=""),
          file=fname_ndx,
          append=TRUE)
  }

  onError <- function() {unlink(c(fname_ndx,
                                  fname_ca_dists,
                                  fname_ca_dists_info,
                                  "dist.xvg"))}
  # compute distances with gromacs
  run.cmd(get.binary("gmx"), c("distance",
                               "-f",
                               traj,
                               "-n",
                               fname_ndx,
                               "-oall",
                               "-select",
                               seq(0, (nrow(ca_pairs)-1))),
          onError=onError)

  # reformat data
  run.cmds(paste("grep -v \"#\" dist.xvg | grep -v \"@\" ",
               "| ",
               get.binary("awk"),
               " '{for(i=2; i <=NF; ++i)",
               "{printf(\" %s\", $i)} printf(\"\\n\")}'",
               " > ",
               fname_ca_dists,
               sep=""),
           onError=onError)

  # remove intermediate file
  unlink("dist.xvg")

  write.caDistances.info(ref     = ref,
                         traj    = traj,
                         mindist = residue.mindist,
                         maxdist = residue.maxdist,
                         fname   = fname_ca_dists_info)

}


#' Generate reaction coordinates.
#'
#' Generate reaction coordinates by selecting columns of the coordinate file(s).
#' An additional <output>.info file is generated giving additional information
#' on which columns were selected.
#'
#' @param coords Character, either single filename or list of filenames of the
#'  coordinate file(s).
#' @param columns Numeric, either single vector or list of vectors indicating
#'  the columns to be selected per coordinate file.
#' @param output Character, name of the output file.
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output files
#'  already exists (default: \code{FALSE}).
#' @export
generate.reactionCoordinates <- function(coords, columns, output, ignoreCache=FALSE) {

  if (xor(is.list(coords), is.list(columns))) {
    stop("Either coordinate files and column selections are both lists, or none is")
  } else if (is.list(coords)) {
    if (length(coords) != length(columns)) {
      stop("List of coordinate files and list of column selections must be of same length.")
    }
  } else {
    coords <- list(coords)
    columns <- list(columns)
  }

  fname_coords <- output
  fname_info   <- paste(fname_coords, ".info", sep="")

  # info file
  info <- do.call(c, lapply(1:length(coords),
                            function(i) {
                              paste(coords[[i]],
                                    paste(columns[[i]], collapse=" "))
                            }))

  # caching
  if (!ignoreCache && file.exists(fname_info) && file.exists(fname_coords)) {
      stored_info <- scan(fname_info, what=character(), sep="\n", quiet=TRUE)

    if (identical(stored_info, info)) {
      return(warning(msg("caching", arg="generate.reactionCoordinates")))
    }
  }

  # write info file
  cat(info, file=fname_info, sep="\n")

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
  cmds <- paste(paste(get.binary("bash"), " -c \"", sep=""),
               get.binary("paste"),
               "-d ' '",
               paste(pipes, collapse=" "),
               ">",
               fname_coords,
               "\" 2>& 1")

  output <- run.cmds(cmds, onError=function() {unlink(c(fname_info, fname_coords))})

  # check for stderr output (due to command substitution the error can not be
  # caught by run.cmds)
  if (length(output) > 0) {
    unlink(c(fname_info, fname_coords))
    stop(paste(output, msg("noExec", arg=cmds), sep="\n"), call.=F)
  }
}

#' Filter a data set
#'
#' Filter a dataset either by selecting every \eqn{N}-th frame or by selecting
#' all frames corresponding to state \eqn{N}.
#'
#' @param coords Character, name of the coordinates file.
#' @param filter Numeric, some integer \eqn{N}.
#'  If \code{states} is \code{NULL} (default), every \eqn{N}-th frame will be selected.
#'  If \code{states} is a filename, N is selected state.
#' @param states Character, name of state trajectory file
#'  (default: \code{NULL}, i.e. not used)
#' @param output Character, name of the ouput file.
#'  If \code{NULL} (default) name is constructed from the input parameters.
#' @param ignoreCache Logical, if \code{TRUE} recompute even if output files
#'  already exists (default: \code{FALSE}).
#' @export
generate.filteredCoordinates <- function(coords, filter, states=NULL, output=NULL, ignoreCache=FALSE) {

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

  if (!ignoreCache && file.exists(output)) {
    warning(msg("caching", arg="generate.filteredCoordinates"))
  } else {
    run.cmds(cmd)
  }
}


