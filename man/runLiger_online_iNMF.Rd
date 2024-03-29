% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCell_Integration.R
\name{runLiger_online_iNMF}
\alias{runLiger_online_iNMF}
\title{Run runLiger_online_iNMF integration (see also http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/online_iNMF_tutorial.html)}
\usage{
runLiger_online_iNMF(
  object,
  liger.k = 30,
  liger.miniBatch_size = 5000,
  liger.max.epochs = 5
)
}
\arguments{
\item{object}{The SingCellaR object.}

\item{liger.k}{Inner dimension of factorization–number of metagenes (default 30). A value in the range 20-50 works well for most analyses.}

\item{liger.miniBatch_size}{Total number of cells in each minibatch (default 5000). This is a reasonable default, but a smaller value such as 1000 may be necessary for analyzing very small datasets. In general, minibatch size should be no larger than the number of cells in the smallest dataset.}

\item{liger.max.epochs}{Maximum number of epochs (complete passes through the data). (default 5)}
}
\description{
Run runLiger_online_iNMF integration (see also http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/online_iNMF_tutorial.html)
}
