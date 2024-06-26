% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbcti_package_function.R
\name{calculate_pseudotime}
\alias{calculate_pseudotime}
\title{calculate pseudotime for each cell}
\usage{
calculate_pseudotime(object, start_state_name)
}
\arguments{
\item{object}{a cell_trajectory object}

\item{start_state_name}{index for the state that starts}
}
\value{
values in the pseudotime slot of cell_trajectory object
}
\description{
calculate pseudotime for each cell
}
\examples{
calculate_pseudotime(sample_data, start_state_name = c('1','2'))
}
