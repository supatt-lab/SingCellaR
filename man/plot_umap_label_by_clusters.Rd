% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellPlots.R
\name{plot_umap_label_by_clusters}
\alias{plot_umap_label_by_clusters}
\title{Plot UMAP with identified clusters}
\usage{
plot_umap_label_by_clusters(
  object,
  show_method = c("walktrap", "louvain", "kmeans", "merged_walktrap", "merged_louvain",
    "merged_kmeans"),
  point.size = 2,
  mark.clusters = TRUE,
  mark.font.size = 10,
  mark.font.color = "yellow",
  add.cluster.in.cells = FALSE,
  cluster.font.size = 1,
  cluster.font.color = "yellow"
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{show_method}{The clustering method names used to generate clusters.}

\item{point.size}{The point size. Default 2}

\item{mark.clusters}{is logical. If TRUE, the number of identified clusters will be displayed.}

\item{mark.font.size}{The font size of the cluster name.}

\item{mark.font.color}{The font color of the cluster name.}

\item{add.cluster.in.cells}{is logical. If TRUE, the cluster name will be displayed in each data point.}

\item{cluster.font.size}{The font size of cluster name in each data point.}

\item{cluster.font.color}{The font color of cluster name in each data point.}
}
\description{
Plot UMAP with identified clusters
}
