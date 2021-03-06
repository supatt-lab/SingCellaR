% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellRNASeq.R
\name{get_variable_genes_by_fitting_GLM_model}
\alias{get_variable_genes_by_fitting_GLM_model}
\title{Identified variable genes by fitting the GLM model using mean vs coefficient of variation}
\usage{
get_variable_genes_by_fitting_GLM_model(
  object,
  mean_expr_cutoff = 0.1,
  disp_zscore_cutoff = 0.1,
  quantile_genes_expr_for_fitting = 0.2,
  quantile_genes_cv2_for_fitting = 0.8
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{mean_expr_cutoff}{The mean expression cutoff.}

\item{disp_zscore_cutoff}{The dispersion of z-score cutoff.}

\item{quantile_genes_expr_for_fitting}{The quantile gene expression used for fitting the model.}

\item{quantile_genes_cv2_for_fitting}{The quantile coefficient of variation used for fitting the model.}
}
\description{
Identified variable genes by fitting the GLM model using mean vs coefficient of variation
}
