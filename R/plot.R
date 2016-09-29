
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
#' @import ggplot2 Biobase edgeR stats
#'
#' @examples
plotExpression <- function(x, cluster = TRUE, scale = FALSE, cpm = FALSE, prior.count = 2, log = TRUE) {
  if (class(x) == "DGEList")
    x <- edgeR::cpm(x, prior.count = prior.count, log = log)

  if (class(x) == "EList")
    x <- x$E

  if (inherits(x, "eSet"))
    x <- exprs(x)

  # Is this needed? I had two implementations, the one with this didn't have the call to cpm above for DGEList... Rethink.
  if (cpm)
    x <- edgeR::cpm(x, prior.count = prior.count, log = log)

  if (cluster) {
    h <- hclust(as.dist(1 - cor(t(x))))
    x <- x[h$order,]
  }

  if (scale)
    x <- t(scale(t(x)))

  d <- reshape2::melt(x, varnames = c("protein", "sample"))

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
  x <- as.matrix(x)
  mode(x) <- "integer"
  x <- x[order(x[,1], x[,2], x[,3]),]
  x <- reshape2::melt(x, varnames = c("protein", "coeficient"))
  x[,3] <- as.factor(x[,3])
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
  d <- reshape2::melt(x$p.value, varnames = c("protein", "coefficient"))
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

  d <- reshape2::melt(y, varnames = c("protein", "sample"))
  d$group <- group[d$sample]
  ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25) + facet_wrap(~group) + theme(aspect.ratio = 1)
}



#' Plot jittered points by group.
#'
#' @param x an object from which a matrix can be obtained.
#' @param selection character or numeric vector with a selection to plot.
#' @param groupCol grouping variable.
#' @param cpm whether to compute cpm().
#' @param prior.count argument for cpm()
#' @param log argument for cpm()
#'
#' @return
#' @export
#'
#' @examples
plotPoints <- function(x, selection = NULL, groupCol = "group", cpm = FALSE, prior.count = 2, log = TRUE) {
  if (cpm) {
    y <- edgeR::cpm(x, prior.count = prior.count, log = log)
  } else {
    if (class(x) == "DGEList")
      y <- x$counts
    if (class(x) == "EList")
      y <- x$E
    if (class(x) == "ExpressionSet")
      y <- exprs(x)
  }

  if (class(x) == "DGEList") {
    group <- x$samples[[groupCol]]
    names(group) <- rownames(x$samples)
  }

  if (class(x) == "EList") {
    group <- x$targets[[groupCol]]
    names(group) <- rownames(x$targets)
  }

  if (class(x) == "ExpressionSet") {
    s <- sampleNames(x)
    group <- x[[groupCol]]
    names(group) <- sampleNames(x)
  }

  if (!is.null(selection))
    y <- y[selection, , drop = FALSE]

  d <- reshape2::melt(y, varnames = c("protein", "sample"))
  d$group <- group[d$sample]

  dd <- d %>% group_by_("protein", "group") %>% summarize_(mean = "mean(value)")

  d %>% ggplot(aes_string(x = "group", y = "value", color = "group")) +
    geom_jitter(width = .3, size = 1, height = 0) +
    facet_wrap(~protein) +
    guides(color = guide_legend(title = "group")) +
    geom_hline(aes_string(yintercept = "mean", color = "group"), data = dd)
}


#' Plot a heatmap with the significance of terms from enrichment results.
#'
#' @param ... results from enrichment, passed down to getEnrichmentResults()
#' @param cutoff p.value cutoff used to filter results.
#' @param what what to plot: P.Up (default) or P.Down. Passed down to getEnrichmentResults()
#' @param use.name logical; whether to return the term's name or the ids.  Passed down to getEnrichmentResults()
#'
#' @return
#' @export
#'
#' @examples
plotEnrichment <- function(..., cutoff = 0.05, what = "P.Up", use.name = TRUE) {
  k <- getEnrichmentResults(..., what = what, use.name = use.name)
  k <- k[rowSums(k < cutoff) > 0,]

  h <- hclust(as.dist(1 - cor(t(k))))
  k <- k[h$order, ]

  d <- k %>% reshape2::melt(varnames = c("term", "group"), value.name = "p.value")

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
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
plotVenn <- function(x, ...) {
  m2l <- function(x) {
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
    l
  }
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid::grid.newpage()
  grid::grid.draw(VennDiagram::venn.diagram(m2l(x), filename = NULL, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", ...))
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

  d <- reshape2::melt(m, varnames = c("sample_i", "sample_j"), value.name = "correlation")

  ggplot(d, aes_string(x = "sample_i", y = "sample_j", fill = "correlation")) +
    geom_tile() +
    labs(x = "", y = "", title = title) +
    viridis::scale_fill_viridis(limits = c(0, 1), guide = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 1, axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}
