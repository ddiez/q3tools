#' Create a matrix with several objects resulting from enrichment analysis with kegga() and goana()
#'
#' @param ... results from enrichment.
#' @param what what to plot: P.Up (default) or P.Down.
#'
#' @return
#' @export
#'
#' @examples
getKEGGResults <- function(..., what = "P.Up") {
  l <- list(...)
  x <- matrix(NA, ncol = length(l), nrow = nrow(l[[1]]))
  rownames(x) <- rownames(l[[1]])
  colnames(x) <- names(l)
  for(i in seq_len(length(l))) {
    x[, i] <- l[[i]][, what]
  }
  x
}

#' Create a matrix with several objects resulting from enrichment analysis with kegga() and goana()
#'
#' @param ... results from enrichment
#' @param what what to plot: P.Up (default) or P.Down.
#' @param use.name logical; whether to return the term's name or the ids.
#'
#' @return
#' @export
#'
#' @examples
getEnrichmentResults <- function(..., what = "P.Up", use.name = FALSE) {
  l <- list(...)
  x <- matrix(NA, ncol = length(l), nrow = nrow(l[[1]]))
  if (use.name)
    rownames(x) <- l[[1]][,1]
  else
    rownames(x) <- rownames(l[[1]])
  colnames(x) <- names(l)
  for(i in seq_len(length(l))) {
    x[, i] <- l[[i]][, what]
  }
  x
}

#' Group-wise missing value imputation
#'
#' @description By default MNAR missing values are assigned 0. MAR values are imputed with some of the methods.
#'
#' @param x matrix with missing values to impute.
#' @param group grouping (column) variable.
#' @param do.mar logical; whether to perform MAR imputation.
#' @param do.mnar logical; whether to perform MNAR imputation.
#' @param mnar.default numerical; default value assigned to MNAR missing values.
#' @param method method used for imputation.
#' @param ... additional arguments passed to the imputation methods.
#'
#' @return
#' @export
#'
#' @examples
imputeGroup <- function(x, group = NULL, do.mar = TRUE, do.mnar = TRUE, mnar.default = 0, method = "mle", ...) {
  s <- split(x, rep(group, each = nrow(x)))
  l <- lapply(s, function(z) {
    tmp <- matrix(z, nrow = nrow(x), byrow = FALSE)

    sel.mnar <- rowSums(is.na(tmp)) == ncol(tmp)

    if (do.mnar)
      tmp[sel.mnar,] <- 0

    if (do.mar) {
      switch(method,
             mle = {
               # MLE
               s <- norm::prelim.norm(tmp[!sel.mnar,])
               th <- norm::em.norm(s, ...)
               tmp[!sel.mnar,] <- norm::imp.norm(s, th, tmp[!sel.mnar,])
             },
             nkk = {
               # KNN
               tmp[!sel.mnar,] <- impute::impute.knn(tmp[!sel.mnar,], ...)$data
             }
      )
    }
    tmp
  })

  m <- do.call(cbind, l)
  rownames(m) <- rownames(x)
  colnames(m) <- colnames(x)
  m
}


#' Get plotting data from plotMDS() in the limma package
#'
#' @param ... pass function arguments to plotMDS.
#'
#' @return
#' @export
#'
#' @importFrom grDevices dev.off png
#'
#' @examples
getMDS <- function(...){
  ff <- tempfile()
  png(filename=ff)
  res <- limma::plotMDS(...)
  dev.off()
  unlink(ff)
  res
}
