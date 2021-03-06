% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellPlots.R
\name{plot_tsne_label_by_AUCell_score}
\alias{plot_tsne_label_by_AUCell_score}
\title{Plot TSNE with AUCell score}
\usage{
plot_tsne_label_by_AUCell_score(
  object,
  geneSets.AUC,
  show_gene_sets = c(),
  isUseCutOffScore = F,
  cutOffScore = 0,
  isSubtract.random.gene = T,
  point.size = 1,
  point.color1 = "gray",
  point.color2 = "red"
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{geneSets.AUC}{The dataframe of AUCell score per gene set.}

\item{show_gene_sets}{The vector of gene set names}

\item{isUseCutOffScore}{is logical. If TRUE, the cutoff score will be applied.}

\item{cutOffScore}{The cutoff score. Default 0}

\item{isSubtract.random.gene}{is logical. If TRUE, the AUC score will be subtracted by the AUCell score from the random gene set.}

\item{point.size}{The point size. Default 1}

\item{point.color1}{The color for a low AUCell score. Default gray}

\item{point.color2}{The color for a high AUCell score . Default red}
}
\description{
Plot TSNE with AUCell score
}
