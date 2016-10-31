#' plotPairs
#'
#' Plots a matrix of pairwise scatter plots (similar to graphics::pairs and GGally::ggpairs).
#'
#' @param x matrix to plot.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotPairs <- function(x) {
  n <- ncol(x)

  comp_matrix <- as.matrix(expand.grid(seq_len(n), seq_len(n)))

  glist <- lapply(seq_len(nrow(comp_matrix)), function(k) {
    z <- comp_matrix[k, ]
    type <- levels(factor(sign(z[1] - z[2])))

    switch(type,
       "0" = {g <- ggplotGrob(qplot(x[,z[1]], geom = "histogram"))},
      "-1" = {g <- ggplotGrob(qplot(x[,z[1]], x[, z[2]]))},
       "1" = {g <- ggplotGrob(qplot(x[,z[1]], x[, z[2]], col = I("red")))}
    )

    gtable_filter(g, "panel")
  })

  gmatrix <- matrix(glist, nrow = n, ncol = n)

  gt <- gtable_matrix("ggpairs", gmatrix, unit(rep(1, n), "null"), unit(rep(1, n), "null"))

  gt <- gtable_add_rows(gt, unit(.2, "null"), pos = 0)
  for (j in 1:ncol(x)) {
    gt <- gtable_add_grob(gt, textGrob(colnames(x)[j], x = unit(.5, "npc"), y = unit(.5, "npc")), t = 1, l = j)
  }

  gt <- gtable_add_cols(gt, unit(.2, "null"), pos = 0)
  for (j in 1:ncol(x)) {
    gt <- gtable_add_grob(gt, textGrob(colnames(x)[j], x = unit(.5, "npc"), y = unit(.5, "npc"), rot = 90), t = j + 1, l = 1)
  }
  gt <- gtable_add_row_space(gt, unit(.1, "null"))
  gt <- gtable_add_col_space(gt, unit(.1, "null"))
  gt <- gtable_add_padding(gt, unit(.1, "null"))


  grid::grid.newpage()
  grid::grid.draw(gt)
  invisible(gt)
}
