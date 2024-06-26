% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbcti_package_function.R
\name{feature_selection_knn}
\alias{feature_selection_knn}
\title{Select features via NMF step 2}
\usage{
feature_selection_knn(object, k = 10, feature)
}
\arguments{
\item{object}{a cell_trajectory object}

\item{k}{number of nearest neighbour}

\item{feature}{the specifeid features}
}
\value{
values in the selected_feature slot of cell_trajectory object
}
\description{
Select features via NMF step 2
}
\examples{
feature_selection_knn(sample_data,feature = c('1','2'), k = 5)
}
