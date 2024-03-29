% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellRNASeq.R
\name{identifyDifferentialGenes_for_all_clusters}
\alias{identifyDifferentialGenes_for_all_clusters}
\title{Identify differentially expressed genes for all identified clusters}
\usage{
identifyDifferentialGenes_for_all_clusters(
  object,
  method = c("wilcoxon"),
  cluster.type = c("louvain", "walktrap", "kmeans", "merged_walktrap", "merged_louvain",
    "merged_kmeans"),
  min.log2FC = 0.2,
  min.expFraction = 0.2
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{method}{The differential gene expression testing method. Default 'wilcoxon'  test}

\item{cluster.type}{The type of clustering method.}

\item{min.log2FC}{The minimum log2FC cutoff. Default 0.20}

\item{min.expFraction}{The minimum fraction of expressing cell frequency cutoff. Default 0.2}
}
\description{
Identify differentially expressed genes for all identified clusters
}
