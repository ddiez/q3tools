#' A bunch of shabby tools for the omics data analyst
#'
#' These package contains some plotting tools at the moment. They help make my
#' life easier when analyzing omics data. But they are not neccesarily the best
#' way to do things. Indeed, I expect them to contain many errors, bad choices
#' and limitations. You have been warned.
#'
#' @name q3tools
#'
#' @import Gviz ggplot2 dplyr tidyr tibble gtable grid
#' @importFrom edgeR cpm
#' @importFrom Biobase exprs pData fData sampleNames
#' @importFrom stats as.dist cor hclust
#' @importFrom reshape2 melt
#' @importFrom IRanges IRanges disjointBins start end width
#' @importFrom grDevices dev.off png
#' @importFrom futile.logger flog.threshold
#'
NULL