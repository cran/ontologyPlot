% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attributes.R
\name{colour_by_population_frequency}
\alias{colour_by_population_frequency}
\title{Function to assign colours to terms based on population frequency of terms}
\usage{
colour_by_population_frequency(
  ontology,
  terms,
  frequencies,
  colour_palette = colorRampPalette(c("Yellow", "Green", "#0099FF"))(10),
  max_colour_freq = max(terms_freq),
  min_colour_freq = min(terms_freq)
)
}
\arguments{
\item{ontology}{\code{ontology_index} object}

\item{terms}{Character vector of ontological terms}

\item{frequencies}{Numeric vector of term frequencies named by term IDs}

\item{colour_palette}{Character vector of colours for the different information contents of the terms to be plotted, going from rare to common}

\item{max_colour_freq}{Numeric value in [0, 1] giving the maximum frequency (to which the dullest color will be assigned)}

\item{min_colour_freq}{Numeric value in [0, 1] giving the minimum frequency (to which the brightest color will be assigned)}
}
\value{
Character vector of colours, named by term
}
\description{
Function to assign colours to terms based on population frequency of terms
}
\seealso{
\code{\link{colour_by_term_set}}, \code{\link{colour_by_frequency}}
}
