% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbcti_package_function.R
\name{distribution_estimation}
\alias{distribution_estimation}
\title{Estimate distribution}
\usage{
distribution_estimation(object, ndraw = 1000, expansion = 1.5, ...)
}
\arguments{
\item{object}{a cell_trajectory object}

\item{ndraw}{number of draws}

\item{expansion}{eapansion parameter for fitted distribution}

\item{...}{cluster index}
}
\value{
values in the distribution_estimation slot of cell_trajectory object
}
\description{
Estimate distribution
}
\examples{
distribution_estimation(sample_data, ndraw = 50, expansion = 1.5, ... = 1,2,3)
}
