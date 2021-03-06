% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellRNASeq.R
\name{add_cell_cycle_genes_score}
\alias{add_cell_cycle_genes_score}
\title{Add the cell cycle gene scores into the cell metadata}
\usage{
add_cell_cycle_genes_score(
  object,
  gmt.file = c(),
  gene_sets = c("G2M_Core", "S_phase_Core")
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{gmt.file}{The GMT file that contains singature gene sets.}

\item{gene_sets}{The vector of gene set names.}
}
\description{
Add the cell cycle gene scores into the cell metadata
}
