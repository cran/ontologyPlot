\name{ontologyPlot}
\docType{package}
\title{
	Functions for Visualising Sets of Ontological Terms
}
\description{
Functions for visualising sets of ontological terms using the `graphviz' layout system.
}
\details{
\tabular{ll}{
Package: \tab ontologyPlot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-01-11\cr
License: \tab GPL (>= 2)\cr
}
This package succeeds the package \code{hpoPlot} with an improved interface and focusing on general ontologies. The key function is \code{\link{onto_plot}}, which creates an object of class \code{ontology_plot} which can be displayed as a graph or exported to \code{dot} format.
}
\author{
Daniel Greene <dg333@cam.ac.uk>

Maintainer: Daniel Greene <dg333@cam.ac.uk>
}
\references{
	Greene D, Richardson S, Turro E (2017). `ontologyX: a suite of R packages for working with ontological data. _Bioinformatics_, 33(7), 1104--1106.

	`The Human Phenotype Ontology project: linking molecular biology and disease through phenotype data', Nucl. Acids Res. (1 January 2014) 42 (D1): D966-D974 doi:10.1093/nar/gkt1026
	Westbury, S. K. et al. (2015). Human Phenotype Ontology annotation and cluster analysis to un- ravel genetic defects in 707 cases with unexplained bleeding and platelet disorders. Genome Medicine. 7 (2015)
}
\keyword{ontology}
\keyword{HPO}
\keyword{MPO}
\keyword{GO}
