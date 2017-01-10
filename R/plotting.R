
#' plot Ramachandran plot
#' @param resno Residue number.
#' @export
plt.ramachandran <- function(resno) {
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
plt.pcaOverview <- function(pca, pcs, corr=FALSE) {
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
#'
#' Takes all available PCAs and plots their cumulative fluctuations in a single plot.
#' @export
plt.cumFlucts <- function() {
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
plt.autocorr <- function(coords, columns, corrlength=0.25) {
  #TODO: finish
}
