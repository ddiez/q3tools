#' Plots a heatmap of a matrix
#'
#' Plots a heatmap with row and/or col clustering if specifief.
#'
#' @param x object to plot.
#' @param row.cluster logical specifying whether to apply hierarchichal clustering to the rows.
#' @param col.cluster logical specifying whether to apply hierarchichal clustering to the columns.
#' @param scale logical specifying whether to scale the rows.
#'
#' @return This function returns no value but has the side effect of producing a plot.
#' @export
#'
#' @docType methods
#' @rdname plotHeatmap-methods
#' @examples
#' NULL
setGeneric("plotHeatmap", function(x, row.cluster = TRUE, col.cluster = FALSE, scale = FALSE)
  standardGeneric("plotHeatmap")
)

#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,matrix-method
setMethod("plotHeatmap", "matrix",
function(x, row.cluster = TRUE, col.cluster = FALSE, scale = FALSE) {
  x <- as.matrix(x)

  if (row.cluster) {
    h <- hclust(as.dist(1 - cor(t(x))))
    x <- x[h$order,]
  }

  if (col.cluster) {
    h <- hclust(as.dist(1 - cor(x)))
    x <- x[, h$order]
  }

  if (scale)
    x <- t(scale(t(x)))

  d <- melt(x, varnames = c("gene", "sample"))
  d$gene <- factor(d$gene, levels = rownames(x)) # this may be needed when gene ids are entrezgene (numeric).

  ggplot(d, aes_string(x = "sample", y = "gene", fill = "value")) +
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
})

#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,ExpressionSet-method
setMethod("plotHeatmap", "ExpressionSet",
function(x, row.cluster = TRUE, col.cluster = FALSE, scale = FALSE) {
  plotHeatmap(exprs(x), row.cluster = row.cluster, col.cluster = col.cluster, scale = scale)
})

#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,EList-method
setMethod("plotHeatmap", "EList",
function(x, row.cluster = TRUE, col.cluster = FALSE, scale = FALSE) {
  plotHeatmap(x$E, row.cluster = row.cluster, col.cluster = col.cluster, scale = scale)
})

#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,DGEList-method
setMethod("plotHeatmap", "DGEList",
function(x, row.cluster = TRUE, col.cluster = FALSE, scale = FALSE) {
  plotHeatmap(x$counts, row.cluster = row.cluster, col.cluster = col.cluster, scale = scale)
})
