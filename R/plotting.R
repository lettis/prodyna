

#' \eq{-\log(x)} transformation.
#'
#' Transform \eq{x} to \eq{-\log(x)} for plotting.
#'
#' @param base Numeric, base to which the logarithm is computed.
.reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-",
                           format(base)),
                    trans,
                    inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}


#' Ramachandran plot.
#'
#' Plot the distribution of \eqn{\phi} and \eqn{\psi}} angles of the residue
#' with the given number.
#'
#' @param resno Numeric, residue number.
#' @param dihedrals Character or data.frame, name of the .dih file or data.frame
#'   containing the dihedral angles.
#' @param dihedralsInfo Character, name of the .dih.info file. \cr
#'   This needs to be set only if \code{dihedrals} is the name of a .dih file
#'   and the file name is not given by <dihedrals>.info
#' @param reslabel Character, label of the residue (to be used in the title).
#' @return ggplot object
#' @import ggplot2
#' @export
plt.ramachandran <- function(resno, dihedrals, dihedralsInfo=NULL, reslabel=NULL) {

  if (is.character(dihedrals)) {
    dih <- read.dihedrals(dihedrals, dihedralsInfo, resno)
  } else {
    dih <- dihedrals
  }
  phi <- dih[[paste("phi", resno, sep="")]]
  psi <- dih[[paste("psi", resno, sep="")]]

  p <- ggplot(data.frame(phi, psi)) +
       stat_bin2d(aes(x=phi, y=psi), bins=180) +
       scale_fill_distiller(palette="YlGnBu",
                            trans=.reverselog_trans()) +
       xlim(-180,180) +
       ylim(-180,180) +
       theme_bw() +
       theme(legend.position = "none")

  if (is.null(reslabel)) {
    p <- p + ggtitle(paste("residue", resno))
  } else {
    p <- p + ggtitle(paste("residue", reslabel))
  }
  return(p)
}

#' plot matrix
#'
#' Plots a matrix as colored tiles.
#' @param x Either a filename with the matrix in table format or a matrix object.
#' @param diverge Use diverging color palette.
#' @param fancy Fancy plotting with interactive elements (default: FALSE).
#' @param zlim Limits for color scale, given as 2D vector.
#'             default: NULL, i.e. min/max values of matrix.
#' @export
plt.matrix <- function(x, diverge=FALSE, fancy=FALSE, zlim=NULL) {
  suppressMessages(require(ggplot2))
  suppressMessages(require(dplyr))
  suppressMessages(require(plotly))
  if (is.character(x)) {
    # interpret as filename
    M <- data.matrix(read.table(x))
  } else {
    M <- data.matrix(x)
  }
  colnames(M) <- 1:ncol(M)
  rownames(M) <- 1:nrow(M)

  if (diverge) {
    clr_palette <- "RdYlBu"
  } else {
    clr_palette <- "YlGnBu"
  }
  # plot rows along y, columns along x
  p <- ggplot(reshape2::melt(M)) +
          geom_raster(aes(y=Var1, x=Var2, fill=value)) +
          scale_y_reverse(breaks=1:nrow(M), expand=c(0,0)) +
          scale_x_continuous(breaks=1:ncol(M), expand=c(0,0)) +
          theme_bw() +
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title = element_blank())
  # set z-axis limits (if given)
  if (is.null(zlim)) {
    p <- p + scale_fill_distiller(palette=clr_palette)
  } else {
    p <- p + scale_fill_distiller(palette=clr_palette,
                                  limits=zlim)
  }
  # fancy plotting for Rnotebooks
  if (fancy) {
    p <- ggplotly(p, tooltip="value") %>%
          config(displayModeBar=FALSE) %>%
          layout(xaxis=list(fixedrange=TRUE)) %>%
          layout(yaxis=list(fixedrange=TRUE))
  }

  p
}

#' PCA overview
#'
#' Plot 2d-proj, 1d-proj and eigenvector content for given PCA.
#'
#' @param fname Character, filename of original coordinates.
#' @param pcs Numerical vector, PC indices.
#' @param corr Logical, if \code{TRUE} use correlation-based PCA (default:
#'   \code{FALSE}).
#' @import ggplot2 GGally data.table
#' @export
plt.pcaOverview <- function(fname, pcs, corr=FALSE) {
  if (corr) {
    proj <- fread(sprintf("%s.projn", fname),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
    vecs <- fread(sprintf("%s.vecn", fname),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
  } else {
    proj <- fread(sprintf("%s.proj", fname),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
    vecs <- fread(sprintf("%s.vec", fname),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
  }
  vec_names <- names(vecs)
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
                                              trans=.reverselog_trans()),
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

#' plot 2d-proj for given PCA
#'
#' plots a 2d free energy landscape (in kT) for the given PCA projection.
#' @param pca Either a PCA reference from the project (e.g. prodyna::projectInfo()$dPCAplus), or a filename.
#' @param dim1 First dimension to plot (default: 1).
#' @param dim2 Second dimension to plot (default: 2).
#' @param corr Plot projections of correlation based PCA. Only works with project reference (default: FALSE).
#' @param diverge Use a diverging color scale to emphasize differences.
#'                Default is FALSE.
#' @export
plt.pcaProj <- function(pca, dim1=1, dim2=2, corr=FALSE, diverge=FALSE) {
  suppressMessages(require(data.table))
  suppressMessages(require(ggplot2))
  suppressMessages(require(dplyr))

  file_read <- function(fname) {
    fread(fname,
          select=c(dim1, dim2),
          verbose=FALSE,
          showProgress=FALSE)
  }

  if (is.data.frame(proj)) {
    proj <- proj %>% select(dim1, dim2)
  } else {
    if (is.list(pca)) {
      .check.projectPath()
      if (corr) {
        .check.filesExist(c(pca$projn, pca$vecn))
        proj <- file_read(pca$projn)
      } else {
        .check.filesExist(c(pca$proj, pca$vec))
        proj <- file_read(pca$proj)
      }
    } else {
      proj <- file_read(pca)
    }
  }

  colnames(proj) <- c("x", "y")

  if (diverge) {
    # diverging color scale
    color_palette <- "RdYlBu"
  } else {
    # linear color scale
    color_palette <- "YlGnBu"
  }

  ggplot(proj) +
    geom_bin2d(bins=200,
               aes(x=x,
                   y=y,
                   fill=-log(..count../max(..count..)))) +
    scale_fill_distiller(palette=color_palette,
                         guide=guide_legend(title="[kT]",
                                            reverse=TRUE)) +
    xlab(paste("PC", dim1)) +
    ylab(paste("PC", dim2)) +
    theme_bw()
}

#' plot cumulative fluctuations
#'
#' Takes all available PCAs and plots their cumulative fluctuations in a single plot.
#' @export
plt.cumFlucts <- function() {
  suppressMessages(require(ggplot2))
  .check.projectPath()
  pd <- projectInfo()

  #TODO different caPCAs!!!

  # check different kinds of PCA for results
  cfs_labels <- c("dPCA+",
                  "dPCA+ (corr)",
                  "caPCA",
                  "caPCA (corr)")
  cfs_filenames <- c(get.fullPath(pd$dPCAplus$val),
                     get.fullPath(pd$dPCAplus$valn),
                     get.fullPath(pd$caPCA$val),
                     get.fullPath(pd$caPCA$valn))
  max_length <- 0
  cfs <- list()
  cfs_labels_selected <- list()
  for (i in 1:length(cfs_labels)) {
    f <- cfs_filenames[i]
    if (file.exists(f) & !dir.exists(f)) {
      cf <- data.table::fread(f, verbose=FALSE, showProgress=FALSE)$V1
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
    scale_x_discrete(limits=1:max_length) +
    #xlim(1, max_length) +
    theme_bw()
}

## plot autocorrelation of observables
##   coords: input file with coordinates
##   columns: vector of column indices to use
##   corrlength: extend of autocorrelation computation.
##               if < 1, ratio of number of data points

#' plot autocorrelation of observables
#'
#' @param coords Filename of coordinates or list of filenames.
#' @param lag.max The maximal lagtime. If < 1 (default: 0.25),
#'                interpret as ratio to total length,
#'                else given in number of timesteps.
#' @param columns Vector with column indices of observables to compute ACF for.
#' @param circular Compute ACF for circular variables (e.g. dihedral angles).
#' @param dt Timestep in [ps]. If given, time axis will be scaled accordingly.
#'           Default: NULL, i.e. express time in number of timesteps.
#' @export
plt.autocor <- function(coords, lag.max=0.25, columns, circular=FALSE, dt=NULL) {
  suppressMessages(require(ggplot2))
  suppressMessages(require(data.table))

  compute_acf <- function(fname) {
    if ( ! file.exists(fname)) {
      fname <- get.fullPath(fname)
    }
    stats.autocor(fname,
                  columns = columns,
                  lag.max = lag.max,
                  circular = circular)
  }
  # compute ACF data (for ensemble: mean and sd)
  if (is.list(coords)) {
    n_files <- length(coords)
    acf_tmp <- lapply(coords, function(fname) {
      do.call("rbind", compute_acf(fname))
    })
    # compute mean
    acf_mean <- acf_tmp[[1]]
    for (i in 2:n_files) {
      acf_mean <- acf_mean + acf_tmp[[i]]
    }
    acf_mean <- acf_mean / n_files
    # compute sigma
    acf_sigma <- 0 * acf_mean
    for (i in 1:n_files) {
      acf_sigma <- (acf_mean-acf_sigma)^2
    }
    acf_sigma <- sqrt(acf_sigma/(n_files-1))
    # reshape and set column names
    acf_mean <- data.frame(t(acf_mean))
    acf_sigma <- data.frame(t(acf_sigma))
    colnames(acf_sigma) <- columns
  } else {
    acf_mean <- compute_acf(coords)
    acf_sigma <- NULL
  }
  # construct time axis
  n <- length(acf_mean[[1]])
  if (is.null(dt)) {
    timeline <-  0:(n-1)
    x_lbl <- "frame"
  } else {
    timeline <- 0:(n-1) * dt
    x_lbl <- "time [ps]"
  }
  # data aggregation
  acf_mean <- data.frame(t=timeline, acf_mean)
  colnames(acf_mean) <- c("t", columns)
  acf_mean <- melt(acf_mean,
                   id.vars=c("t"),
                   value.name="mu")
  if (is.null(acf_sigma)) {
    acf_data <- acf_mean
  } else {
    acf_sigma <- data.frame(t=timeline, acf_sigma)
    colnames(acf_sigma) <- c("t", columns)
    acf_sigma <- melt(acf_sigma,
                      id.vars=c("t"),
                      value.name="sigma")
    acf_data <- dplyr::left_join(acf_mean, acf_sigma, by=c("t", "variable"))
    acf_data <- dplyr::mutate(acf_data,
                              ymin=sapply(mu-sigma, function(x) {max(0.1, x)}),
                              ymax=sapply(mu+sigma, function(x) {max(0.1, x)}))
  }

  if (is.null(acf_sigma)) {
    p <- ggplot(acf_data) +
      geom_line(aes(x=t, y=mu, color=variable))
  } else {
    p <- ggplot(acf_data, aes(t)) +
      geom_ribbon(aes(ymin=ymin,
                      ymax=ymax,
                      fill=variable), alpha=0.2) +
      geom_line(aes(y=mu, color=variable))
  }
  p <- p +
    scale_y_log10(limits=c(0.1, 1)) +
    xlab(x_lbl) +
    ylab("ACF") +
    theme_bw() +
    theme(legend.title = element_blank())

  p
}

#' plot state geometries as Ramacolor plots
#'
#' @param statetraj File with state trajectory
#' @param states State selection for Ramacolor plot. Either a vector with state numbers, or NULL, meaning all (default).
#' @param dihedrals File with dihedral angles. If NULL (default), take dihedrals from project.
#' @export
plt.ramacolor <- function(statetraj, states=NULL, dihedrals=NULL) {
  suppressMessages(require(ggplot2))
  if (is.null(dihedrals)) {
    .check.projectPath()
    pd <- projectInfo()
    dihedrals <- pd$dihedrals
  }
  dih <- data.table::fread(dihedrals,
                           verbose = FALSE,
                           showProgress = FALSE)
  traj <- data.table::fread(statetraj,
                            verbose = FALSE,
                            showProgress = FALSE)[[1]]
  if (is.null(states)) {
    states <- sort(unique(traj))
  }
  n_residues <- ncol(dih) / 2
  # construct Ramacolor dataframe
  rc_data <- do.call("rbind", lapply(states, function(s) {
    sel_dih <- dih[traj == s, ]
    do.call("rbind", lapply(1:n_residues, function(i_res) {
      phis <- sel_dih[[2*i_res-1]]
      psis <- sel_dih[[2*i_res]]
      rgb_vec <- rama2rgb(phis, psis)
      data.frame(ires=i_res, state=s, r=rgb_vec[1], g=rgb_vec[2], b=rgb_vec[3])
    }))
  }))

  ggplot(rc_data) +
    geom_raster(aes(x=state, y=ires, fill=rgb(r, g, b)), hjust=0) +
    scale_fill_identity() +
    xlab("state") +
    ylab("residue") +
    xlim(1, length(states)) +
    theme_bw()
}

#' plot state trajectory comparison
#'
#' Compare state trajectories based on their overlap of state populations.
#' Trajectories must be of same length and must have identical state labels.
#' Colored states are taken from traj1 and separate into black states, taken
#' from traj2.
#'
#' @param traj1 Either a vector encoding first state trajectory or
#'              a filename pointing to the trajectory.
#' @param traj2 Either a vector encoding second state trajectory or
#'              a filename pointing to the trajectory.
#' @export
plt.stateTrajComparison <- function(traj1, traj2) {
  check_traj <- function(traj) {
    if (is.character(traj)) {
      traj <- data.table::fread(traj,
                                verbose=FALSE,
                                showProgress=FALSE)[[1]]
    }
    traj
  }
  traj1 <- check_traj(traj1)
  traj2 <- check_traj(traj2)
  states <- sort(unique(c(traj1, traj2)))
  n_states <- length(states)
  overlap <- matrix(nrow=n_states, ncol=n_states)
  idx <- seq(1, n_states)
  for (i in idx) {
    state_is_i <- traj2[(traj1==states[i])]
    for (j in idx) {
      state_is_ij <- (state_is_i == states[j])
      overlap[i,j] <- sum(state_is_ij)
    }
  }
  rownames(overlap) <- paste("", states)
  colnames(overlap) <- paste(" ", states)
  #TODO: nicer plot, no scaling, etc
  circlize::chordDiagram(overlap,
                         grid.col=c(rainbow(n_states),
                                    rep("black", n_states)))
}



