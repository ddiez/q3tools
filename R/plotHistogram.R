#' Plot a hisogram of values.
#'
#' Plot a histogram, maybe with facets as a function of some grouping variable.
#'
#' @param x     object from which to extract a matrix to plot the histogram.
#' @param group grouping variable (columns).
#'
#' @return
#' @export
#'
#' @docType methods
#' @rdname plotHistogram-methods
#' @examples
setGeneric("plotHistogram", function(x, group = NULL)
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

#' @rdname plotHistogram-methods
#' @aliases plotHistogram,ExpressionSet-method
setMethod("plotHistogram", "ExpressionSet",
function(x, group = NULL) {
  plotHistogram(exprs(x), group = group)
})

#' @rdname plotHistogram-methods
#' @aliases plotHistogram,EList-method
setMethod("plotHistogram", "EList",
  function(x, group = NULL) {
    plotHistogram(x$E, group = group)
})

#' @rdname plotHistogram-methods
#' @aliases plotHistogram,DGEList-method
setMethod("plotHistogram", "DGEList",
function(x, group = NULL) {
  plotHistogram(x$counts, group = group)
})
