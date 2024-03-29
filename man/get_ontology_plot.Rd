% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphs.R
\name{get_ontology_plot}
\alias{get_ontology_plot}
\title{Get \code{ontology_plot} object}
\usage{
get_ontology_plot(
  ontology,
  terms,
  edge_attributes = list(color = "#000000", lty = "solid"),
  ...
)
}
\arguments{
\item{ontology}{\code{ontology_index} object}

\item{terms}{Character vector of ontological terms}

\item{edge_attributes}{List of properties to set for arrows (note, these properties will be used for all arrow).}

\item{...}{Named graphical parameters. These must either be vectors of values the same length as \code{terms}, or of length 1 if they should be used for all terms.}
}
\value{
\code{ontology_plot} object.
}
\description{
Function to create \code{ontology_plot} objects where all graphical parameters to be used must be specified.
}
