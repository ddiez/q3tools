#' Create a matrix with several objects resulting from enrichment analysis with kegga() and goana()
#'
#' @param ... results from enrichment
#' @param what what to plot: P.Up (default) or P.Down.
#' @param ontology for GO results what ontology to use (default BP).
#' @param use.name logical; whether to return the term's name or the ids.
#' @param short.names whether to use abbreviate() to shorten rownames (handy for long GO descriptions).
#' @param minlength minimum length for abbreviate().
#'
#' @return a matrix with p-values from several enrichment analyses.
#' @export
#'
#' @examples
#' NULL
getEnrichmentResults <- function(..., what = "P.Up", ontology = "BP", use.name = FALSE, short.names = FALSE, minlength = 40) {
  l <- list(...)

  # this assumes all the objects passed are homogeneous (i.e. all GO or all KEGG)
  if ("Ont" %in% colnames(l[[1]])) {
    l <- lapply(l, function(ll) {
      ll %>% filter_("Ont == ontology")
    })
  }
  x <- matrix(NA, ncol = length(l), nrow = nrow(l[[1]]))
  if (use.name)
    rownames(x) <- l[[1]][,1]
  else
    rownames(x) <- rownames(l[[1]])
  if (short.names)
    rownames(x) <- abbreviate(rownames(x), minlength = minlength)
  colnames(x) <- names(l)
  for(i in seq_len(length(l))) {
    tmp <- l[[i]]
    x[, i] <- tmp[, what]
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
#' @return a matrix with missing values imputed.
#' @export
#'
#' @examples
#' NULL
imputeGroup <- function(x, group = NULL, do.mar = TRUE, do.mnar = TRUE, mnar.default = 0, method = "knn", ...) {
  if (do.mnar) {
    d <- dimnames(x)
    s <- split(x, rep(group, each = nrow(x)))
    l <- lapply(s, function(z) {
      tmp <- matrix(z, nrow = nrow(x), byrow = FALSE)

      sel.mnar <- rowSums(is.na(tmp)) == ncol(tmp)
      if (any(sel.mnar))
        tmp[sel.mnar, ] <- mnar.default
      tmp
    })
    x <- do.call(cbind, l)
    dimnames(x) <- d
  }

  if (do.mar) {
    switch(method,
           mle = {
             # MLE
             s <- norm::prelim.norm(x)
             th <- norm::em.norm(s, ...)
             x <- norm::imp.norm(s, th, x)
           },
           knn = {
             # KNN
             x <- impute::impute.knn(data = x, ...)$data
           }
    )
  }
  x
}


#' getMDS
#'
#' Get plotting data from plotMDS() in the limma package
#'
#' @param ... pass function arguments to plotMDS.
#'
#' @return a list with information about MDS, including plotting data.
#' @export
#'
#' @examples
#' NULL
getMDS <- function(...){
  limma::plotMDS(..., plot = FALSE)
}
