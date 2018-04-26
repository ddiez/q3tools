#' A bunch of shabby tools for the omics data analyst
#'
#' This package contains some plotting tools at the moment. They help make my
#' life easier when analyzing omics data. But they are not neccesarily the best
#' way to do things. Indeed, I expect them to contain many errors, bad choices
#' and limitations. You have been warned.
#'
#' @name q3tools
#' @docType package
#'
#' @import tibble gtable
#' @rawNamespace import(ggplot2, except = "exprs")
#' @importFrom dplyr %>% filter_ group_by_ summarize_ select
#' @importFrom grid grid.newpage grid.draw unit textGrob
#' @importFrom Gviz IdeogramTrack GenomeAxisTrack BiomartGeneRegionTrack DataTrack plotTracks chromosome
#' @importFrom edgeR cpm
#' @importFrom limma plotMDS
#' @importFrom Biobase exprs pData fData sampleNames
#' @importFrom stats as.dist cor hclust
#' @importFrom reshape2 melt
#' @importFrom IRanges IRanges disjointBins start end width ranges
#' @importFrom grDevices dev.off png
#' @importFrom futile.logger flog.threshold
#'
NULL
