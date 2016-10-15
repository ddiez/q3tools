#' Plot a hisogram of values.
#'
#' Plot a histogram, maybe with facets as a function of some grouping variable.
#'
#' @param x     object from which to extract a matrix to plot the histogram.
#' @param group grouping variable (columns).
#' @param ... further arguments passed to methods.
#'
#' @return
#' @export
#'
#' @docType methods
#' @rdname plotHistogram-methods
#' @examples
setGeneric("plotHistogram", function(x, group = NULL, ...)
  standardGeneric("plotHistogram")
)

#' @rdname plotHistogram-methods
#' @aliases plotHistogram,matrix-method
setMethod("plotHistogram", "matrix",
function(x, group = NULL) {
  d <- melt(x, varnames = c("gene", "sample"))
  if (!is.null(group)) {
    names(group) <- colnames(group)
    d$group <- group[d$sample]
  }
  g <- ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25)

  if (!is.null(group))
    g + facet_wrap(~group)
  else
    g
})

#' @param groupCol Column name from which extract grouping variable.
#'
#' @rdname plotHistogram-methods
#' @aliases plotHistogram,ExpressionSet-method
setMethod("plotHistogram", "ExpressionSet",
function(x, group = NULL, groupCol = NULL) {
  y <- exprs(x)

  if (!is.null(group)) {
    s <- sampleNames(x)
    names(group) <- sampleNames(x)
  } else {
    if (!is.null(groupCol)) {
      group <- x[[groupCol]]
    }
  }
})

# plotHistogram <- function(x) {
#   if (class(x) == "DGEList") {
#     y <- x$counts
#     group <- x$samples$group
#     names(group) <- rownames(x$samples)
#   }
#
#   if (class(x) == "EList") {
#     y <- x$E
#     group <- x$targets$group
#     names(group) <- rownames(x$targets)
#   }
#
#   if (class(x) == "ExpressionSet") {
#     y <- exprs(x)
#     s <- sampleNames(x)
#     group <- x$group
#     names(group) <- sampleNames(x)
#   }
#
#   d <- melt(y, varnames = c("protein", "sample"))
#   d$group <- group[d$sample]
#   ggplot(d, aes_string(x = "value")) + geom_histogram(bins = 25) + facet_wrap(~group) + theme(aspect.ratio = 1)
# }
