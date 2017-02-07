#' Plots the mean-variance trend.
#'
#' This function uses the information stored in the object returned by voom()
#' when called with save.plot = TRUE.
#'
#' @param x Output of voom() function.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotVoom <- function(x) {
  d <- x$voom.xy
  l <- x$voom.line
  qplot(d$x, d$y, size = I(.1), xlab = d$xlab, ylab = d$ylab, main = "voom: Mean-variance trend") + annotate("line", x = l$x, y = l$y, color = "red", size = 1)
}

#' Plot a heatmap of TestResult object from limma package
#'
#' This plots a heatmap of the TestResult object returned by decideTests. The
#' columns are the coefficients and the rows are the genes.
#'
#' @param x TestResult object.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotResult <- function(x) {
  x <- unclass(x) # For TestResult objects. Make this more general.

  x <- x[do.call(order, as.list(as.data.frame(x))),] # to sort by each column sequentially.

  x <- melt(x, varnames = c("protein", "coefficient"))
  x$value <- factor(x$value)

  ggplot(x, aes_string(x="coeficient", y = "protein", fill = "value")) +
    geom_raster() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("blue3", "white", "red3"))
}


#' plotPvalue
#'
#' Plot distribution of p-values by group.
#'
#' @param x object with a p.value element.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotPvalue <- function(x) {
  d <- melt(x$p.pvalue, varnames = c("protein", "coefficient"))
  ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25) + facet_wrap(~coefficient) + labs(title = "Distribution of p-value") + theme(aspect.ratio = 1)
}

#' plotEnrichment
#'
#' Plot a heatmap with the significance of terms from enrichment results.
#'
#' @param ... results from enrichment, passed down to getEnrichmentResults()
#' @param cutoff p.value cutoff used to filter results.
#' @param what what to plot: P.Up (default) or P.Down. Passed down to getEnrichmentResults()
#' @param ontology for GO results what ontology to use (default BP).
#' @param use.name logical; whether to return the term's name or the ids.  Passed down to getEnrichmentResults()
#' @param short.names whether to use abbreviate() to shorten rownames (handy for long GO descriptions).
#' @param minlength minimum length for abbreviate().
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotEnrichment <- function(..., cutoff = 0.05, what = "P.Up", ontology = "BP", use.name = TRUE, short.names = FALSE, minlength = 40) {
  k <- getEnrichmentResults(..., what = what, ontology = ontology, use.name = use.name, short.names = short.names, minlength = minlength)
  k <- k[rowSums(k < cutoff) > 0,]

  h <- hclust(as.dist(1 - cor(t(k))))
  k <- k[h$order, ]

  d <- melt(k, varnames = c("term", "group"), value.name = "p.value")

  ggplot(d, aes_string(x = "group", y = "term", fill = "-log10(p.value)")) +
    geom_raster() + viridis::scale_fill_viridis(guide = "legend") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = what) +
    theme(axis.ticks.y = element_blank())
}


#' plotVenn
#'
#' Plot a venn diagram with VennDiagram package using some nice defaults
#'
#' @param x matrix of values (similar to what is fed to VennDiagram in limma package).
#' @param ... additional parameters passed down to venn.diagrama()
#' @param add.universe logical; whether to add a group containing all elements.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
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
  grid::grid.draw(
    VennDiagram::venn.diagram(
      m2l(x, add.universe = add.universe),
      filename = NULL,
      fontfamily = "sans",
      cat.fontfamily = "sans",
      main.fontfamily = "sans",
      ...
    )
  )
}


#' plotMds
#'
#' ggplo2-style MDS plot using plotMDS() in the limma package to get the data.
#'
#' @param x an object from which a data matrix can be obtained.
#' @param group grouping variable.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
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


#' plotCorrelation
#'
#' Plot a heatmap of a correlation matrix.
#'
#' Computes the correlation matrix and passes it to plotHeatmap() with
#' appropriate arguments.
#'
#' @param x and object from which an matrix can be obtained.
#' @param title title of the plot.
#' @param cluster logical; whether to cluster rows/columns.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotCorrelation <- function(x, title = "Sample correlation", cluster = FALSE) {
  suppressMessages(
    plotHeatmap(cor(x), row.cluster = cluster, col.cluster = cluster, scale = FALSE) +
      theme(axis.text.y = element_text()) +
      labs(x = "sample_i", y = "sample_j", title = title) +
      guides(fill = guide_legend("correlation", reverse = TRUE)) +
      viridis::scale_fill_viridis(limit = c(0, 1))
  )
}

#' Basic gene ploting function
#'
#' Plots a gene from a specified assembly using Gviz package. Optionaly plots some data (e.g. read alignments). Enables zooming.
#'
#' @param symbol symbol of the gene to plot. Will be used to search using BiomartGeneRegionTrack.
#' @param genome genome assembly (e.g. mm10 or hg38).
#' @param add.ideogram whether to add an IdeogramTrack (FALSE by default to spead up testing).
#' @param add.data list with data to be added as a DataTrack.
#' @param from start genomic coordinates.
#' @param to end genomic coordinates.
#' @param biomart mart object (optional).
#' @param ... further arguments passed to plotTracks.
#'
#' @return nothing but produces a plot as side effect.
#' @export
#'
#' @examples
#' NULL
plotGene <- function(symbol, genome, add.ideogram = FALSE, add.data = NULL, from = NULL, to = NULL, biomart = NULL, ...) {
  itrack <- NULL
  if (add.ideogram)
    itrack <- IdeogramTrack(genome = genome)
  atrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE)
  if (!is.null(biomart)) {
    bmtrack <- BiomartGeneRegionTrack(symbol = symbol,
                                      biomart = biomart,
                                      name = "EnsEMBL",
                                      geneSymbols = TRUE,
                                      showTranscriptID = TRUE,
                                      col.line = NULL,
                                      col = NULL,
                                      col.title = "black",
                                      background.title = "white")
  } else {
    bmtrack <- BiomartGeneRegionTrack(symbol = symbol,
                                      genome = genome,
                                      name = "EnsEMBL",
                                      geneSymbols = TRUE,
                                      showTranscriptID = TRUE,
                                      col.line = NULL,
                                      col = NULL,
                                      col.title = "black",
                                      background.title = "white")
  }


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

    if (!is.null(from) || !is.null(to)) {
      if (is.null(from))
        from <- min(start(bmtrack))

      if (is.null(to))
        to <- max(end(bmtrack))
    }
  }

  tl <- c(
    itrack,
    atrack,
    bmtrack,
    dtrack
  )
  plotTracks(tl, chromosome = chromosome(bmtrack), from = from, to = to, ...)
  invisible(bmtrack)
}

