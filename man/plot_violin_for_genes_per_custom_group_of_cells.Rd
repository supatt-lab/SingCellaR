% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellPlots.R
\name{plot_violin_for_genes_per_custom_group_of_cells}
\alias{plot_violin_for_genes_per_custom_group_of_cells}
\title{Plot violin-plot for gene expression per custom group of cells}
\usage{
plot_violin_for_genes_per_custom_group_of_cells(
  object,
  custom_group_of_cells = list(),
  gene_list = c(),
  take_log2 = T,
  xlab.text.size = 5,
  point.size = 0.2,
  point.alpha = 0.1,
  grid.ncol = 3,
  grid.nrow = 3
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{custom_group_of_cells}{The list of custome group of cells.}

\item{gene_list}{The vector of gene names.}

\item{take_log2}{is logical. If TRUE, log2 expression will be applied.}

\item{xlab.text.size}{The font size of the label on the x-axis. Default 5}

\item{point.size}{The size of data point. Default 0.2}

\item{point.alpha}{The alpha parameter of the data point. Default 0.1}

\item{grid.ncol}{the column number of the grid of the containing plots to be displayed. Default 3}

\item{grid.nrow}{the row number of the grid of the containing plots to be displayed. Default 3}
}
\description{
Plot violin-plot for gene expression per custom group of cells
}
