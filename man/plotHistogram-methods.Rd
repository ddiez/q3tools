% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHistogram.R
\docType{methods}
\name{plotHistogram}
\alias{plotHistogram}
\alias{plotHistogram,matrix-method}
\alias{plotHistogram,ExpressionSet-method}
\alias{plotHistogram,EList-method}
\alias{plotHistogram,DGEList-method}
\title{Plot a hisogram of values.}
\usage{
plotHistogram(x, group = NULL)

\S4method{plotHistogram}{matrix}(x, group = NULL)

\S4method{plotHistogram}{ExpressionSet}(x, group = NULL)

\S4method{plotHistogram}{EList}(x, group = NULL)

\S4method{plotHistogram}{DGEList}(x, group = NULL)
}
\arguments{
\item{x}{object from which to extract a matrix to plot the histogram.}

\item{group}{grouping variable (columns).}
}
\description{
Plot a histogram, maybe with facets as a function of some grouping variable.
}
\examples{
NULL
}
