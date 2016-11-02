#' Plot ranges avoiding overlapping
#'
#' This function plots a set of ranges in their locations. Overlapping is
#' avoided. This function is based on a similar one in the IRanges vignette,
#' but using ggplot2 graphics instead.
#'
#' @param x an IRanges object or some object that can be cohereced to IRanges.
#' @param sep vertical separation between ranges.
#' @param ... further arguments passed down to methods.
#'
#' @return NULL
#' @export
#' @rdname plotRanges-methods
#'
#' @examples
#' NULL
setGeneric("plotRanges", function(x, sep = 0.5, ...)
  standardGeneric("plotRanges")
)

#' @rdname plotRanges-methods
#' @aliases plotRanges,IRanges-method
#' @importFrom IRanges IRanges disjointBins start end width
setMethod("plotRanges", "IRanges",
function(x, sep = 0.5) {
  # determine y-distance.
  height <- 1
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  ybottom <- bins * (sep + height) - height
  d <- x %>% as.data.frame %>% select(-width)
  colnames(d) <- c("xmin", "xmax")
  d <- cbind(d, ymin = ybottom, ymax = ybottom + height) %>% rownames_to_column("range")

  # compute plot limits and set breaks.
  r <- range(x)
  b <- seq(start(r), end(r))

  ggplot(d,
         aes_string(
           xmin = "xmin",
           xmax = "xmax",
           ymin = "ymin",
           ymax = "ymax",
           fill = "range",
           color = "range",
           label = "range"
         )) +
    geom_vline(xintercept = b, lty = "dotted", color = "grey", size = .5) +
    geom_rect() +
    geom_text(aes_string(x = "xmin + (xmax - xmin) / 2", y = "ymin + (ymax - ymin) / 2"), color = "black") +
    scale_x_continuous(breaks = b) +
    labs(x = "position", y = "") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
})

#' @rdname plotRanges-methods
#' @aliases plotRanges,GRanges-method
#' @param seqnames sequence name of the ranges.
setMethod("plotRanges", "GRanges",
function(x, sep = 0.5, seqnames = NULL) {
  if (is.null(seqnames)) stop("seqnames must not be NULL.")
  plotRanges(ranges(x[seqnames(x) == seqnames]), sep = sep)
})
