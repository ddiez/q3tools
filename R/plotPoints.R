#' Plot jittered points by group from a matrix.
#'
#' Takes a matrix and produces a dots plot using geom jitter.
#'
#' @details Alternatively points are grouped by the groups assigned to the
#' group argument, which is applied to the columns. Each row is plotted in a
#' different panel using facet_wrap. The label of the facets is by default
#' the rownames but it can be passed using the argument label. The scales
#' of the plots are by default the same but this can be controled with the
#' scale argument (e.g. changed to 'free_y').
#'
#' @param x a matrix or an object that can be coherce to a matrix.
#' @param group grouping variable (columns).
#' @param label custom labels (rows).
#' @param scales scales argument passed down to facet_wrap (default: fixed).
#'
#' @return NULL
#' @export
#'
#' @docType methods
#' @rdname plotPoints-methods
#'
#' @examples
#' NULL
setGeneric("plotPoints", function(x, group = NULL, label = NULL, scales = "fixed")
  standardGeneric("plotPoints")
)

#' @rdname plotPoints-methods
#' @aliases plotPoints,matrix-method
setMethod("plotPoints", "matrix",
function(x, group = NULL, label = NULL, scales = "fixed") {
  x <- as.matrix(x)

  if (is.null(colnames(x)))
    colnames(x) <- seq_len(1:ncol(x))

  if (is.null(rownames(x)))
    rownames(x) <- seq_len(1:nrow(x))

  if (!is.null(group))
    names(group) <- colnames(x) # colnames can be NULL!

  if (is.null(label)) {
    label <- rownames(x)
    names(label) <- rownames(x)
  } else {
    names(label) <- rownames(x)
  }

  d <- melt(x, varnames = c("gene", "sample"))
  d$gene <- factor(d$gene, levels = rownames(x)) # this may be needed when gene ids are entrezgene (numeric).
  d$group <- group[as.character(d$sample)]

  dd <- d %>% group_by_("gene", "group") %>%
    summarize_(mean = "mean(value, na.rm = TRUE)")

  d %>% ggplot(aes_string(x = "group", y = "value", color = "group")) +
    geom_jitter(width = .5, size = 1, height = 0) +
    facet_wrap(~gene, scales = scales, labeller = as_labeller(label)) +
    guides(color = guide_legend(title = "group")) +
    geom_hline(aes_string(yintercept = "mean", color = "group"), data = dd)
})

#' @rdname plotPoints-methods
#' @aliases plotPoints,matrix-method
setMethod("plotPoints", "ExpressionSet",
function(x, group = NULL, label = NULL, scales = "fixed") {
  plotPoints(exprs(x), group = group, label = label, scales = scales)
})

#' @rdname plotPoints-methods
#' @aliases plotPoints,EList-method
setMethod("plotPoints", "EList",
function(x, group = NULL, label = NULL, scales = "fixed") {
  plotPoints(x$E, group = group, label = label, scales = scales)
})


#' @rdname plotPoints-methods
#' @aliases plotPoints,DGEList-method
setMethod("plotPoints", "DGEList",
function(x, group = NULL, label = NULL, scales = "fixed") {
  plotPoints(x$counts, group = group, label = label, scales = scales)
})
