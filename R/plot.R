
#' Plots a heatmap of a matrix
#'
#' @param x object to plot.
#' @param cluster logical specifying whether to apply hierarchichal clustering to the rows.
#' @param scale logical specifying whether to scale the rows.
#' @param cpm logical specifying
#' @param prior.count numeric value passed to cpm().
#' @param log logical passed to cpm() to specified whether values should be log-transformed.
#'
#' @return This function returns no value but has the side effect of producing a plot.
#' @export
#' @import ggplot2 dplyr tidyr tibble
#' @importFrom edgeR cpm
#' @importFrom Biobase exprs pData fData sampleNames
#' @importFrom stats as.dist cor hclust
#'
#' @examples
plotExpression <- function(x, cluster = TRUE, scale = FALSE, cpm = FALSE, prior.count = 2, log = TRUE) {
  if (class(x) == "DGEList")
    x <- cpm(x, prior.count = prior.count, log = log)

  if (class(x) == "EList")
    x <- x$E

  if (inherits(x, "eSet"))
    x <- exprs(x)

  # Is this needed? I had two implementations, the one with this didn't have the call to cpm above for DGEList... Rethink.
  if (cpm)
    x <- cpm(x, prior.count = prior.count, log = log)

  if (cluster) {
    h <- hclust(as.dist(1 - cor(t(x))))
    x <- x[h$order,]
  }

  if (scale)
    x <- t(scale(t(x)))

  d <- melt(x, rows.name = "protein", cols.name = "sample")

  ggplot(d, aes_string(x = "sample", y = "protein", fill = "value")) +
    geom_raster() +
    labs(x = "sample", y = "protein", title = "Protein expression") +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = .5),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    ) +
    viridis::scale_fill_viridis(guide = guide_legend(reverse = TRUE)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
}

#' Plots the mean-variance trend as computed by voom().
#'
#' @param x Output of voom() function (with save.plot = TRUE)
#'
#' @return
#' @export
#'
#' @examples
plotVoom <- function(x) {
  d <- x$voom.xy
  l <- x$voom.line
  qplot(d$x, d$y, size = I(.1), xlab = d$xlab, ylab = d$ylab, main = "voom: Mean-variance trend") + annotate("line", x = l$x, y = l$y, color = "red", size = 1)
}

#' Title
#'
#' @param x TestResult object.
#'
#' @return
#' @export
#'
#' @examples
plotResult <- function(x) {
  x <- unclass(x) # For TestResult objects. Make this more general.

  x <- x[do.call(order, as.list(as.data.frame(x))),] # to sort by each column sequentially.

  x <- melt(x, rows.name = "protein", cols.name = "coeficient")
  x$value <- factor(x$value)

  ggplot(x, aes_string(x="coeficient", y = "protein", fill = "value")) +
    geom_raster() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("blue3", "white", "red3"))
}


#' Title
#'
#' @param x object with a p.value element.
#'
#' @return
#' @export
#'
#' @examples
plotPvalue <- function(x) {
  d <- melt(x$p.value, rows.name = "protein", cols.name = "coefficient")
  ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25) + facet_wrap(~coefficient) + labs(title = "Distribution of p-value") + theme(aspect.ratio = 1)
}


#' Title
#'
#' @param x object from which to extract a matrix to plot the histogram.
#'
#' @return
#' @export
#'
#' @examples
plotHistogram <- function(x) {
  if (class(x) == "DGEList") {
    y <- x$counts
    group <- x$samples$group
    names(group) <- rownames(x$samples)
  }

  if (class(x) == "EList") {
    y <- x$E
    group <- x$targets$group
    names(group) <- rownames(x$targets)
  }

  if (class(x) == "ExpressionSet") {
    y <- exprs(x)
    s <- sampleNames(x)
    group <- x$group
    names(group) <- sampleNames(x)
  }

  d <- melt(y, rows.name = "protein", cols.name = "sample")
  d$group <- group[d$sample]
  ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25) + facet_wrap(~group) + theme(aspect.ratio = 1)
}



#' Plot jittered points by group.
#'
#' @param x an object from which a matrix can be obtained.
#' @param selection character or numeric vector with a selection to plot.
#' @param group character vector with grouping variable.
#' @param groupCol grouping variable name (to obtain group from x).
#' @param cpm whether to compute cpm().
#' @param prior.count argument for cpm().
#' @param log argument for cpm().
#' @param scales argument passed to facet_wrap().
#'
#' @return
#' @export
#'
#' @examples
plotPoints <- function(x, selection = NULL, group = NULL, groupCol = "group", cpm = FALSE, prior.count = 2, log = TRUE, scales = "fixed") {
  if (cpm) {
    y <- edgeR::cpm(x, prior.count = prior.count, log = log)
  } else {
    if (class(x) == "DGEList")
      y <- x$counts
    if (class(x) == "EList")
      y <- x$E
    if (class(x) == "ExpressionSet")
      y <- exprs(x)
    if (class(x) == "matrix")
      y <- x
  }

  if (class(x) == "DGEList" && is.null(group)) {
    group <- x$samples[[groupCol]]
    names(group) <- rownames(x$samples)
  }

  if (class(x) == "EList" && is.null(group)) {
    group <- x$targets[[groupCol]]
    names(group) <- rownames(x$targets)
  }

  if (class(x) == "ExpressionSet" && is.null(group)) {
    s <- sampleNames(x)
    group <- x[[groupCol]]
    names(group) <- sampleNames(x)
  }

  if (class(x) == "matrix") {
    if (is.null(group)) stop("group must not be NULL for matrix objects.")
    names(group) <- colnames(x)
  }

  if (!is.null(selection))
    y <- y[selection, , drop = FALSE]

  d <- melt(y, rows.name = "protein", cols.name = "sample")
  d$group <- group[d$sample]

  dd <- d %>% group_by_("protein", "group") %>% summarize_(mean = "mean(value, na.rm = TRUE)")

  d %>% ggplot(aes_string(x = "group", y = "value", color = "group")) +
    geom_jitter(width = .5, size = 1, height = 0) +
    facet_wrap(~protein, scales = scales) +
    guides(color = guide_legend(title = "group")) +
    geom_hline(aes_string(yintercept = "mean", color = "group"), data = dd)
}


#' Plot a heatmap with the significance of terms from enrichment results.
#'
#' @param ... results from enrichment, passed down to getEnrichmentResults()
#' @param cutoff p.value cutoff used to filter results.
#' @param what what to plot: P.Up (default) or P.Down. Passed down to getEnrichmentResults()
#' @param ontology for GO results what ontology to use (default BP).
#' @param use.name logical; whether to return the term's name or the ids.  Passed down to getEnrichmentResults()
#'
#' @return
#' @export
#'
#' @examples
plotEnrichment <- function(..., cutoff = 0.05, what = "P.Up", ontology = "BP", use.name = TRUE) {
  k <- getEnrichmentResults(..., what = what, ontology = ontology, use.name = use.name)
  k <- k[rowSums(k < cutoff) > 0,]

  h <- hclust(as.dist(1 - cor(t(k))))
  k <- k[h$order, ]

  d <- k %>% melt(rows.name = "term", cols.name = "group", value.name = "p.value")

  ggplot(d, aes_string(x = "group", y = "term", fill = "-log10(p.value)")) +
    geom_raster() + viridis::scale_fill_viridis(guide = "legend") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = what) +
    theme(axis.ticks.y = element_blank())
}


#' Plot a venn diagram with VennDiagram package using some nice defaults
#'
#' @param x matrix of values (similar to what is fed to VennDiagram in limma package).
#' @param ... additional parameters passed down to venn.diagrama()
#' @param add.universe logical; whether to add a group containing all elements.
#'
#' @return
#' @export
#'
#' @importFrom futile.logger flog.threshold
#'
#' @examples
plotVenn <- function(x, add.universe = FALSE, ...) {
  m2l <- function(x, add.universe = FALSE) {
    l <- lapply(seq_len(ncol(x)), function(i) {
      if (is.null(rownames(x)))
        (1:nrow(x))[x[,i]]
      else
        rownames(x)[x[,i]]
    })
    if (is.null(colnames(x)))
      names(l) <- paste0("group-", 1:length(l))
    else
      names(l) <- colnames(x)
    if (add.universe)
      l$total <- rownames(x)
    l
  }
  flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid::grid.newpage()
  grid::grid.draw(VennDiagram::venn.diagram(m2l(x, add.universe = add.universe), filename = NULL, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", ...))
}


#' ggplo2-style MDS plot using plotMDS() in the limma package to get the data.
#'
#' @param x an object from which a data matrix can be obtained.
#' @param group grouping variable.
#'
#' @return
#' @export
#'
#' @examples
plotMds <- function(x, group = NULL) {
  if(class(x) == "RangedSummarizedExperiment")
    x <- SummarizedExperiment::assay(x)

  x <- getMDS(x)

  d <- data.frame(x = x$x, y = x$y, sample = names(x$x))
  d$group <- group
  d$index <- as.character(1:nrow(d))

  g <- ggplot(d, aes_string(x = "x", y = "y")) +
    theme(aspect.ratio = 1) +
    labs(x = paste(x$axislabel, 1), y = paste(x$axislabel, 2))

  if (is.null(group))
    g + geom_text(aes_string(label = "index"), size = 3)
  else
    g + geom_text(aes_string(label = "index", color = "group"), size = 3)
}


#' Plot a heatmap of a correlation matrix.
#'
#' @param x and object from which an matrix can be obtained.
#' @param title title of the plot.
#' @param cluster logical; whether to cluster rows/columns.
#'
#' @return
#' @export
#'
#' @examples
plotCorrelation <- function(x, title = "Sample correlation", cluster = FALSE) {
  if (class(x) == "ExpressionSet")
    x <- exprs(x)

  if(class(x) == "RangedSummarizedExperiment")
    x <- SummarizedExperiment::assay(x)

  m <- cor(x)

  if (cluster) {
    h <- hclust(as.dist(1 - cor(m)))
    m <- m[h$order, h$order]
  }

  d <- melt(m, rows.name = "sample_i", cols.name = "sample_j", value.name = "correlation")

  ggplot(d, aes_string(x = "sample_i", y = "sample_j", fill = "correlation")) +
    geom_tile() +
    labs(x = "", y = "", title = title) +
    viridis::scale_fill_viridis(limits = c(0, 1), guide = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 1, axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}

#' Basic gene ploting function
#'
#' @param symbol symbol of the gene to plot. Will be used to search using BiomartGeneRegionTrack
#' @param genome genome assembly (e.g. mm10 or hg38)
#' @param add.ideogram whether to add an IdeogramTrack (FALSE by default to spead up testing)
#' @param ... further arguments passed to plotTracks
#'
#' @return nothing but produces a plot as side effect.
#' @export
#'
#' @import Gviz
#'
#' @examples
plotGene <- function(symbol, genome, add.ideogram = FALSE, add.data = NULL, ...) {
  itrack <- NULL
  if (add.ideogram)
    itrack <- IdeogramTrack(genome = genome)
  atrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE)
  bmtrack <- BiomartGeneRegionTrack(symbol = symbol,
                                    genome = genome,
                                    name = "EnsEMBL",
                                    geneSymbols = TRUE,
                                    showTranscriptID = TRUE,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title = "black",
                                    background.title = "white")

  dtrack <- NULL
  if (!is.null(add.data)) {

    if (!is.list(add.data))
      add.data <- list(add.data)

    if (is.null(names(add.data)))
      names(add.data) <- paste0("DataTrack-", seq_len(length(add.data)))

    dtrack <- lapply(names(add.data), function(n) {
      x <- add.data[[n]]
      DataTrack(x,
                name = n,
                type = "histogram",
                col.histogram = "cornflowerblue",
                fill.histogram = "cornflowerblue",
                col.title = "black",
                background.title = "white")
    })
  }

  tl <- c(
    itrack,
    atrack,
    bmtrack,
    dtrack
  )
  plotTracks(tl, chromosome = chromosome(bmtrack), ...)
}

