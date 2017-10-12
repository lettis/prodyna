#' Negative Logarithm.
#'
#' Transform \eqn{x} to \eqn{-log(x)} for plotting.
#'
#' @param base Numeric, base to which the logarithm is computed
#'  (default: \eqn{e})
#' @importFrom scales trans_new log_breaks
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv   <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)),
            trans,
            inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

#' Ramachandran plot.
#'
#' Plot the distribution of \eqn{\phi} and \eqn{\psi} angles for a given residue.
#'
#' @param resno Numeric, residue number.
#' @param dihedrals Character or data.frame, name of the .dih file or data.frame
#'   containing the dihedral angles.
#' @param reslabel Character, label of the residue (to be used in the title).
#' @return ggplot object
#' @import ggplot2
#' @export
plt.ramachandran <- function(resno, dihedrals, reslabel=NULL) {

  if (is.character(dihedrals)) {
    dih <- read.dihedrals(dihedrals     = dihedrals,
                          resnos        = resno,
                          dihedralsInfo = NULL)
  } else {
    dih <- dihedrals
  }

  if (is.null(reslabel)) {
    title <- paste("residue", resno)
  } else {
    title <- paste("residue", reslabel)
  }

  phi <- dih[[paste("phi", resno, sep="")]]
  psi <- dih[[paste("psi", resno, sep="")]]

  p <- ggplot(data.frame(phi, psi)) +
       stat_bin2d(aes(x=phi, y=psi), bins=180) +
       scale_fill_distiller(palette="YlGnBu", trans=reverselog_trans()) +
       xlim(-180,180) +
       ylim(-180,180) +
       theme_bw() +
       theme(legend.position = "none") +
       xlab(expression(phi)) +
       ylab(expression(psi)) +
       ggtitle(title)

  return(p)
}

#' Matrix plot.
#'
#' Plots a matrix as colored tiles.
#'
#' @param x Character or matrix object, filename of table-formatted file or a
#'  matrix object.
#' @param diverge Logical, should a diverging color palette be used?
#'   (default:\code{FALSE}).
#' @param fancy Logical, should fancy plotting with interactive elements be
#'   performed? (default: \code{FALSE}).
#' @param zlim Numeric vector, minimum and maximum value for color scale.
#'   If \code{NULL} (default) limits are given by the min/max values of the matrix.
#' @import ggplot2
#' @importFrom plotly ggplotly config layout
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 melt
#' @export
plt.matrix <- function(x, diverge=FALSE, fancy=FALSE, zlim=NULL) {
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
  p <- ggplot(melt(M)) +
          geom_raster(aes(y=Var1, x=Var2, fill=value)) +
          scale_y_reverse(breaks=1:nrow(M), expand=c(0,0)) +
          scale_x_continuous(breaks=1:ncol(M), expand=c(0,0)) +
          theme_bw() +
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title=element_blank())
  # set z-axis limits (if given)
  if (is.null(zlim)) {
    p <- p + scale_fill_distiller(palette=clr_palette)
  } else {
    p <- p + scale_fill_distiller(palette=clr_palette, limits=zlim)
  }
  # fancy plotting for Rnotebooks
  if (fancy) {
    p <- ggplotly(p, tooltip="value") %>%
          config(displayModeBar=FALSE) %>%
          layout(xaxis=list(fixedrange=TRUE)) %>%
          layout(yaxis=list(fixedrange=TRUE))
  }
  return(p)
}

#' PCA overview
#'
#' Plot 2d-proj, 1d-proj and eigenvector content for given PCA.
#'
#' @param coords Character, name of coordinates file for which PCA has been
#'  performed.
#' @param pcs Numeric vector, PC indices.
#' @param corr Logical, if \code{TRUE} use correlation-based PCA (default:
#'   \code{FALSE}).
#' @import ggplot2
#' @importFrom data.table fread
#' @importFrom GGally ggpairs putPlot
#' @export
plt.pcaOverview <- function(coords, pcs, corr=FALSE) {
  if (corr) {
    fname_proj <- sprintf("%s.projn", coords)
    fname_vecs <- sprintf("%s.vecn", coords)
  } else {
    fname_proj <- sprintf("%s.proj", coords)
    fname_vecs <- sprintf("%s.vec", coords)
  }

  proj <- fread(fname_proj, select=pcs, verbose=FALSE, showProgress=FALSE)
  vecs <- fread(fname_vecs, select=pcs, verbose=FALSE, showProgress=FALSE)

  vec_names <- names(vecs)
  n_dih     <- dim(vecs)[1]
  vecs <- do.call("rbind",
                  lapply(vec_names, function(n) {
                    data.frame(dih=vecs[[n]], v=n, ndx=1:n_dih)
                  }))
  plt <- ggpairs(proj, lower="blank", upper="blank")
  seq <- 1:ncol(proj)

  for (x in seq){
    for (y in seq) {

      if (x < y) {
        # lower
        subplt <- ggplot(proj, aes_string(x=names(proj)[x],
                                          y=names(proj)[y])) +
                  stat_bin2d(bins=200) +
                  scale_fill_distiller(palette="YlGnBu",
                                       trans=reverselog_trans())
      } else if (y < x) {
        # upper
        subplt <- ggplot(vecs[vecs$v==vec_names[x] | vecs$v==vec_names[y],],
                      aes(x=ndx, y=abs(dih), group=v, color=v)) +
                  geom_line(size=2)
      } else {
        # diag
        subplt <- ggplot(proj, aes_string(x=names(proj)[x])) +
                  stat_bin(bins=200,
                           aes(y=-log(..density..)),
                           geom="line",
                           position="identity")
      }
    plt <- putPlot(plt, subplt, y, x)
    }
  }
  return(plt)
}

#' Plot a 2d-projection for a given PCA.
#'
#' plots a 2d free energy landscape (in kT) for the given PCA projection.
#' @param projected Character or data.frame, filename of the .proj file or data
#'  frame containing the projected data.
#' @param dim1 Numeric, index of PC to be used as first dimension, i.e. column
#'  to select from file or data frame (default: 1).
#' @param dim2 Numeric, index of PC to be used as second dimension, i.e. column
#'  to select from file or data frame (default: 2).
#' @param diverge Logical, should a diverging color scale be used?
#'  (default: \code{FALSE})
#' @import ggplot2
#' @importFrom data.table fread
#' @export
plt.pcaProj <- function(projected, dim1=1, dim2=2, diverge=FALSE) {

  if (is.character(projected)) {
    pcs <- fread(projected, select=c(dim1, dim2), verbose=F, showProgress=F)
  } else {
    pcs <- projected[, c(dim1, dim2)]
  }

  colnames(pcs) <- c("x", "y")

  if (diverge) {
    # diverging color scale
    color_palette <- "RdYlBu"
  } else {
    # linear color scale
    color_palette <- "YlGnBu"
  }

  ggplot(pcs) +
    geom_bin2d(bins=200,
               aes(x=x,
                   y=y,
                   fill=-log(..count../max(..count..)))) +
    scale_fill_distiller(palette = color_palette,
                         guide   = guide_colourbar(title="[kT]",
                                            reverse=TRUE)) +
    xlab(paste("PC", dim1)) +
    ylab(paste("PC", dim2)) +
    theme_bw()
}

#' Plot cumulative fluctuations.
#'
#' Takes all available PCAs and plots their cumulative fluctuations in a single plot.
#' @param eigenvals Character vector, names of .val/.valn files holding the
#'  eigenvalues.
#' @param labels Character vector, PCA method descriptions corresponding to the
#'  respective files in \code{eigenvals}.
#' @import ggplot2
#' @importFrom data.table fread
#' @export
plt.cumFlucts <- function(eigenvals, labels) {

  cfs        <- list()

  if(is.character(eigenvals) && !all(file.exists(eigenvals))) {
    stop(msg("missingFile", eigenvals))
  }

  for (i in 1:length(labels)) {

    if(is.character(eigenvals)){
      cf <- fread(eigenvals[i], verbose=FALSE, showProgress=FALSE)$V1
    } else {
      cf <- eigenvals[[i]]
    }
    cf <- cumsum(cf/sum(cf))

    cfs[[i]]   <- data.frame(1:length(cf), cf, labels[[i]])
    colnames(cfs[[i]]) <- c("PC", "cumfluct", "method")
  }

  df <- do.call("rbind", cfs)

  ggplot(df) +
    geom_line(aes(x=PC, y=cumfluct, color=method), size=1) +
    theme_bw() +
    ylab("cumulative fluctuation")
}

#' Plot autocorrelation.
#'
#' For a single coordinates file the autocorrelation as a function of time lag
#' is plotted for each column.
#' For multiple coordinate files the mean autocorrelation and variance is
#' plotted (for each column the mean/variance is compute over all files).
#'
#' @param coords Character or (list of) data.frame(s), vector of coordinate
#' filenames or a (list of) data frame(s).
#'  If multiple filenames/data frames are given, mean autocorrelation and
#'  variance are plotted.
#' @param lag.max Numeric, maximal lagtime. If < 1 (default: 0.25), this is
#'  interpreted as a ratio to the total length. Otherwise the lagtime is
#'  assumed to be given in number of timesteps/frames.
#' @param columns Numeric vector, column indices of observables to compute
#'  ACF for.
#' @param circular Logical, should ACF be computed for circular variables,
#'  e.g. dihedral angles? (default: \code{FALSE})
#' @param dt Numeric, timestep in [ps]. If not \code{NULL}, the time axis
#'  will be scaled accordingly. Otherwise (default) time is expressed in the
#'  number of timesteps.
#' @param logy Logical, plot with logarithmic y-scale.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join mutate
#' @export
plt.autocor <- function(coords, columns=NULL, lag.max=0.25, circular=F, dt=NULL, logy=F) {

  compute_acf <- function(coord) {
    stats.autocor(coord,
                  columns = columns,
                  lag.max = lag.max,
                  circular = circular)
  }
  # compute ACF data (for ensemble: mean and sd)
  if (is.data.frame(coords)) {
    n_coords <- 1
  } else {
    n_coords <- length(coords)   # might be 1 as well
  }

  if (n_coords > 1) {

    acf_tmp <- lapply(coords, compute_acf)
    # compute mean
    acf_mean <- Reduce('+', acf_tmp)
    acf_mean <- acf_mean/n_coords
    # compute sigma
    acf_sigma <- Reduce('+', lapply(acf_tmp, function(acf) {(acf-acf_mean)^2}))
    acf_sigma <- sqrt(acf_sigma/(n_coords-1))

  } else {
    acf_mean  <- compute_acf(coords)
    acf_sigma <- NULL
  }
  # construct time axis
  n <- length(acf_mean[[1]])
  if (is.null(dt)) {
    timeline <-  0:(n-1)
    x_lbl    <- "time lag in frames"
  } else {
    timeline <- 0:(n-1) * dt
    x_lbl    <- "time lag in ps"
  }

  if(is.null(columns)) {
    columns <- 1:ncol(acf_mean)
  }
  # data aggregation
  acf_mean <- data.frame(t=timeline, acf_mean)
  names(acf_mean) <- c("t", columns)

  acf_mean <- melt(acf_mean, id.vars=c("t"), value.name="mu")

  if (is.null(acf_sigma)) {
    acf_data <- acf_mean
  } else {
    acf_sigma           <- data.frame(t=timeline, acf_sigma)
    colnames(acf_sigma) <- c("t", columns)

    acf_sigma <- melt(acf_sigma, id.vars=c("t"), value.name="sigma")

    acf_data <- left_join(acf_mean, acf_sigma, by=c("t", "variable"))
    acf_data <- mutate(acf_data,
                       ymin = sapply(mu-sigma, function(x) {x}),
                       ymax = sapply(mu+sigma, function(x) {x}))
  }

  if (is.null(acf_sigma)) {
    p <- ggplot(acf_data) + geom_line(aes(x=t, y=mu, color=variable))
  } else {
    p <- ggplot(acf_data, aes(t)) +
      geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=variable), alpha=0.2) +
      geom_line(aes(y=mu, color=variable))
  }
  p <- p +  xlab(x_lbl) +
            ylab("ACF") +
            theme_bw() +
            theme(legend.title = element_blank())


  if (logy) {
    p <- p + scale_y_log10(limits=c(0.01, 1))
  }

  return(p)
}

#' Ramacolor plots of state geometries.
#'
#' Color-coded dihedral content per state and residue.
#'
#' TODO: all residues are read from file and indexed 1:n_residues
#'
#' @param statetraj Character or data frame, filename of state trajectory file
#'  or data frame holding the state trajectory.#' .
#' @param dihedrals Character or data frame, filename of the .dih file or data
#'  frame holding dihedral angles.
#' @import ggplot2
#' @importFrom data.table fread
#' @export
plt.ramacolor <- function(statetraj, dihedrals) {

  if (is.character(statetraj)) {
    traj <- factor(fread(statetraj, verbose = F, showProgress = F)[[1]])
  } else {
    traj <- statetraj
  }

  if (is.character(dihedrals)) {
    dih  <- fread(dihedrals, verbose = F, showProgress = F)
  } else {
    dih <- dihedrals
  }

  states     <- levels(traj)
  n_residues <- ncol(dih)/2

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
    geom_raster(aes(x=state, y=ires, fill=rgb(r, g, b))) +
    scale_fill_identity() +
    scale_x_discrete(labels = states, breaks = states, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("state") +
    ylab("residue")
}

#' Plot state trajectory comparison.
#'
#' Compare state trajectories based on their overlap of state populations.
#' Trajectories must be of same length and must have identical state labels.
#' Colored states are taken from \code{traj1} and separate into black states,
#' taken from \code{traj2}.
#'
#' @param traj1 Either a vector encoding first state trajectory or
#'              a filename pointing to the trajectory.
#' @param traj2 Either a vector encoding second state trajectory or
#'              a filename pointing to the trajectory.
#' @importFrom data.table fread
#' @importFrom circlize chordDiagram circos.par circos.clear
#' @export
plt.stateTrajComparison <- function(traj1, traj2) {
  check_traj <- function(traj) {
    if (is.character(traj)) {
      traj <- fread(traj, verbose=FALSE, showProgress=FALSE)[[1]]
    } else {
      traj
    }
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

  sgap <- 2
  bgap <- 50

  rowsums <- sum(rowSums(abs(overlap)))
  colsums <- sum(colSums(abs(overlap)))

  n_sector <- 2 * n_states

  row_sector_degree <- (360-sgap*(n_sector-2)-bgap*2) * (rowsums/(rowsums+colsums)) +
    sgap*(n_states-1)

  start_degree <- 0 - (180-row_sector_degree)/2

  gaps <- c(rep(sgap, n_states - 1),
           bgap,
           rep(sgap, n_states - 1),
           bgap)

  circos.par(gap.after = gaps, start.degree = start_degree)

  p <- chordDiagram(overlap,
               grid.col=c(rainbow(n_states), rep("black", n_states)),
               annotationTrack =c('name', 'grid'))
  circos.clear()
}


#' Plot sorted per-frame populations for given radii.
#'
#' @param pops Character or data frame, prefix of population files or data frame
#'  with neighbourhood populations per frame where each column  corresponds to a
#'  particular neighbourhood radius.
#' @param radii Numeric vector, selection of radii. If \code{NULL} (default),
#' read all available population files.
#' @param logy Logical, plot with logarithmic y-scale.
#' @param select Numeric, plot only every n-th data point, i.e. for select=100,
#'  the population of every 100-th frame is plotted.
#'  A large value will improve performance significantly.
#'  If \code{NULL} (default) select is chosen such that approx 1000 data points
#'  are plotted.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plt.populations <- function(pops, radii=NULL, logy=TRUE, select=NULL) {

  if(is.character(pops)) {
    pops <- read.populations(pops, radii)
  }

  if(is.null(select)) {
    select <- floor(nrow(pops)/1000)  # plot approx 1000 data points
  }

  for(i in colnames(pops)) {
    pops[[i]] <- sort(pops[[i]], decreasing=TRUE)
  }

  n_frames <- nrow(pops)
  pops$id <- 1:n_frames

  pops <- pops[1:floor(n_frames/select) * select,]

  p <- ggplot(melt(pops, id.vars="id")) +
    geom_line(aes(x=id, y=value, color=variable)) +
    xlab("frames") +
    ylab("sorted populations") +
    scale_color_discrete(name="Radius")

  if (logy) {
    p <- p + scale_y_log10()
  }
  return(p)
}


#' Plot hierarchical network of MPP lumping
#'
#' @param dirname Directory of MPP run
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate count arrange
#' @importFrom igraph graph plot.igraph layout_with_kk
#'
#' @export
plt.mppNetwork <- function(dirname) {
  #### get the data
  # helper to get different data files from MPP directory
  get_mpp_data <- function(fname) {
    data.table::fread(paste(dirname,
                            fname,
                            sep="/"),
                      verbose=FALSE,
                      showProgress=FALSE)
  }
  # get populations at different qmin levels
  mpp_pops <- lapply(list.files(dirname, pattern="mpp_pop*"), function(fname){
    pops <- get_mpp_data(fname)
    colnames(pops) <- c("state", "pop")
    pops
  })
  # helper function to get max pops from union of p1 and p2
  get_max_pops <- function(p1, p2) {
    all_states <- unique(c(p1$state, p2$state))
    pops <- c()
    for (s in all_states) {
      pops <- c(pops, max(p1$pop[p1$state==s],
                          p2$pop[p2$state==s],
                          na.rm=TRUE))
    }
    data.frame(state=all_states, pop=pops)
  }
  # get total max pops for all states (i.e. pops before lumping)
  max_pops <- Reduce(get_max_pops, mpp_pops)
  # get transitions (i.e. at what qmin level are states lumped?)
  transitions <- get_mpp_data("mpp_transitions.dat")
  colnames(transitions) <- c("from", "to", "qmin")
  transitions$qmin <- round(transitions$qmin, digits=4)

  #### construct the network
  pops <- (max_pops %>% arrange(state) %>% select(pop))[[1]]
  g <- graph(edges=as.vector(rbind(transitions$from, transitions$to)))
  g$weights <- 5*transitions$qmin
  plot(g,
       vertex.color="white",
       vertex.size=log(pops),
       edge.arrow.size=0.2,
       edge.width=g$weights,
       layout=layout_with_kk(g))
}

#' Plot waiting time distributions.
#'
#' @param wtd data.frame, waiting time distributions
#' @param max_frame Numeric, limit the x-axis to obtain a nicer plot.
#' @param ploty Logical, plot with logarithmic y-scale.
#' @import ggplot2
#' @export
plt.wtDistributions <- function(wtd, max_frame=NULL, logy=F) {

  # state and window size as factor
  wtd$wsize <- factor(wtd$wsize)
  wtd$state <- factor(wtd$state)

  title <- "Waiting Time Distribution"

  plt <-ggplot(wtd,
               aes(x=frame,
               y=probability,
               group=wsize,
               color=wsize)) +
        geom_line() +
        ggtitle(title) +
        scale_color_discrete(name="window size") +
        facet_grid(state ~ .,
                   labeller = labeller(.rows=function(i) {paste("state", i)}))

  if(!is.null(max_frame)) {
    plt <- plt + xlim(0, max_frame)
  }
  if (logy) {
    plt <- plt + scale_y_log10()
  }
  return(plt)
}

