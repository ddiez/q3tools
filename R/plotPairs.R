#' plotPairs
#'
#' Plots a matrix of pairwise scatter plots a la base::pairs using ggplot2 graphics.
#' Inspired by GGally::ggpairs, but aimed to be simpler.
#'
#' @param x matrix to plot.
#' @param geom.low geom used for lower triangle (default: point)
#' @param geom.mid geom used for diagonal (default: histogram)
#' @param geom.up geom used for upper triangle (default: point)
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotPairs <- function(x, geom.low = "point", geom.mid = "histogram", geom.up = "point") {
  n <- ncol(x)

  comp_matrix <- as.matrix(expand.grid(seq_len(n), seq_len(n)))

  glist <- lapply(seq_len(nrow(comp_matrix)), function(k) {
    z <- comp_matrix[k, ]
    type <- levels(factor(sign(z[1] - z[2])))

    to_function <- function(e) {
      f <- paste0("geom_", e, "()")
      eval(parse(text = f))
    }

    geom_low <- to_function(geom.low)
    geom_mid <- to_function(geom.mid)
    geom_top <- to_function(geom.up)

    switch(type,
       "0" = {g <- ggplotGrob(ggplot(mapping = aes(x[,z[1]])) + geom_mid)},
      "-1" = {g <- ggplotGrob(ggplot(mapping = aes(x[,z[1]], x[, z[2]])) + geom_top)},
       "1" = {g <- ggplotGrob(ggplot(mapping = aes(x[,z[1]], x[, z[2]])) + geom_low)}
    )

    gtable_filter(g, "panel")
  })

  gmatrix <- matrix(glist, nrow = n, ncol = n, byrow = TRUE)

  gt <- gtable_matrix("ggpairs", gmatrix, unit(rep(1, n), "null"), unit(rep(1, n), "null"))

  gt <- gtable_add_rows(gt, unit(1, "lines"), pos = 0)
  for (j in 1:ncol(x)) {
    gt <- gtable_add_grob(gt, textGrob(colnames(x)[j], x = unit(.5, "npc"), y = unit(.5, "npc")), t = 1, l = j)
  }

  gt <- gtable_add_cols(gt, unit(1, "lines"), pos = 0)
  for (j in 1:ncol(x)) {
    gt <- gtable_add_grob(gt, textGrob(colnames(x)[j], x = unit(.5, "npc"), y = unit(.5, "npc"), rot = 90), t = j + 1, l = 1)
  }
  gt <- gtable_add_row_space(gt, unit(.1, "null"))
  gt <- gtable_add_col_space(gt, unit(.1, "null"))
  gt <- gtable_add_padding(gt, unit(.1, "null"))


  grid.newpage()
  grid.draw(gt)
  invisible(gt)
}
