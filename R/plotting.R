

#' transform x -> -log(x) for plotting
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


#' plot Ramachandran plot
#' @param resno Residue number.
#' @param dihedrals File with dihedrals. If NULL (default), choose default
#'                  dihedrals from project management.
#' @export
plt.ramachandran <- function(resno, dihedrals=NULL) {
  if (is.null(dihedrals)) {
    .check.projectPath()
    p <- projectInfo()
    dihedrals <- p$dihedrals
  }
  suppressMessages(require(ggplot2))
  dih <- read.dihedrals(resno, dihedrals)
  phi <- dih[[paste("phi", resno, sep="")]]
  psi <- dih[[paste("psi", resno, sep="")]]
  p <- ggplot(data.frame(phi, psi)) +
    stat_bin2d(aes(x=phi, y=psi), bins=180) +
    scale_fill_distiller(palette="YlGnBu",
                         trans=.reverselog_trans()) +
    xlim(-180,180) +
    ylim(-180,180) +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle(paste("residue", resno))

  p
}

#' plot matrix
#'
#' Plots a matrix as colored tiles.
#' @param x Either a filename with the matrix in table format or a matrix object.
#' @param diverge Use diverging color palette.
#' @param fancy Fancy plotting with interactive elements (default: FALSE).
#' @export
plt.matrix <- function(x, diverge = FALSE, fancy = FALSE) {
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
          scale_y_reverse(breaks=1:nrow(M)) +
          scale_x_continuous(breaks=1:ncol(M)) +
          scale_fill_distiller(palette=clr_palette) +
          theme_bw() +
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank())

  if (fancy) {
    p <- ggplotly(p, tooltip="value") %>%
          config(displayModeBar=FALSE) %>%
          layout(xaxis=list(fixedrange=TRUE)) %>%
          layout(yaxis=list(fixedrange=TRUE))
  }

  p
}


#' plot 2d-proj, 1d-proj and eigenvector content for given PCA
#' @param pca List with filenames, pointing to projections and eigenvectors.
#'            Format (cov): "proj": "coords.proj", "vec": "coords.vec".
#'            Format (cov): "projn": "coords.projn", "vecn": "coords.vecn".
#' @param pcs Vector of PC indices.
#' @param corr Use correlation-based PCA (default: FALSE).
#' @export
plt.pcaOverview <- function(pca, pcs, corr=FALSE) {
  .check.projectPath()
  suppressMessages(require(data.table))
  suppressMessages(require(GGally))
  suppressMessages(require(ggplot2))
  if (corr) {
    proj <- fread(.check.filePath(pca$projn),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
    vecs <- fread(.check.filePath(pca$vecn),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
  } else {
    proj <- fread(.check.filePath(pca$proj),
                  select=pcs,
                  verbose=FALSE,
                  showProgress=FALSE)
    vecs <- fread(.check.filePath(pca$vec),
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

  file_read <- function(fname) {
    fread(fname,
          select=c(dim1, dim2),
          verbose=FALSE,
          showProgress=FALSE)
  }
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

  if (diverge) {
    # diverging color scale
    color_palette <- "RdYlBu"
  } else {
    # linear color scale
    color_palette <- "YlGnBu"
  }

  ggplot(proj) +
    geom_bin2d(bins=200,
               aes(x=V1,
                   y=V2,
                   fill=-log(..count../max(..count..)))) +
    scale_fill_distiller(palette="YlGnBu",
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
#' @param coords Filename of coordinates.
#' @param lag.max The maximal lagtime. If < 1 (default: 0.25), interpret as ratio to total length, else given in number of timesteps.
#' @param columns Vector with column indices of observables to compute ACF for.
#' @param circular Compute ACF for circular variables (e.g. dihedral angles).
#' @param dt Timestep in [ps]. If given, time axis will be scaled accordingly. Default: NULL, i.e. express time in number of timesteps.
#' @export
plt.autocor <- function(coords, lag.max=0.25, columns, circular=FALSE, dt=NULL) {
  suppressMessages(require(ggplot2))
  suppressMessages(require(data.table))
  if ( ! file.exists(coords)) {
    coords <- get.fullPath(coords)
  }
  acf_data <- stats.autocor(coords,
                            columns = columns,
                            lag.max = lag.max,
                            circular = circular)
  n <- length(acf_data[[1]])
  if (is.null(dt)) {
    timeline <-  0:(n-1)
    x_lbl <- "frame"
  } else {
    timeline <- 0:(n-1) * dt
    x_lbl <- "time [ps]"
  }
  acf_data <- data.frame(t=timeline, acf_data)
  colnames(acf_data) <- c("t", columns)
  acf_data <- melt(acf_data, id.vars=c("t"))

  ggplot(acf_data) +
    geom_line(aes(x=t, y=value, color=variable)) +
    scale_y_log10(limits=c(0.1, 1)) +
    xlab(x_lbl) +
    ylab("ACF") +
    theme_bw() +
    theme(legend.title = element_blank())
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
