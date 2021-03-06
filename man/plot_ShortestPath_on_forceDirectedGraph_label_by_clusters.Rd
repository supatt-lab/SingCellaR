% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellPlots.R
\name{plot_ShortestPath_on_forceDirectedGraph_label_by_clusters}
\alias{plot_ShortestPath_on_forceDirectedGraph_label_by_clusters}
\title{Plot the shortest-path on the force-directed graph with identified clusters}
\usage{
plot_ShortestPath_on_forceDirectedGraph_label_by_clusters(
  object,
  start_cluster = "",
  end_cluster = "",
  show_method = c("walktrap", "louvain", "kmeans"),
  vertex.size = 2,
  path.size = 2,
  path.type = 1,
  path.color = "red",
  mark.clusters = TRUE,
  mark.font.size = 10,
  mark.font.color = "yellow",
  add.cluster.in.cells = FALSE,
  cluster.font.size = 1,
  cluster.font.color = "yellow",
  edge.size = 0.2,
  edge.color = "gray"
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{start_cluster}{The starting cluster.}

\item{end_cluster}{The ending cluster.}

\item{show_method}{The clustering method names used to generate clusters.}

\item{vertex.size}{The node size. Default 2}

\item{path.size}{The path size. Default 2}

\item{path.type}{The path type. Default 1}

\item{path.color}{The path color. Default red}

\item{mark.clusters}{is logical. If TRUE, the number of identified clusters will be displayed.}

\item{mark.font.size}{The font size of the cluster name. Default 10}

\item{mark.font.color}{The font color of the cluster name. Default yellow}

\item{add.cluster.in.cells}{is logical. If TRUE, the cluster name will be displayed in each data point.}

\item{cluster.font.size}{The font size of cluster name in each data point. Default 1}

\item{cluster.font.color}{The font color of cluster name in each data point. Default yellow}

\item{edge.size}{The size of the edge. Default 0.2}

\item{edge.color}{The edge color. Default gray}
}
\description{
Plot the shortest-path on the force-directed graph with identified clusters
}
