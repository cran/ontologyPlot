#' Get set of HPO terms appropriate for showing in a grid
#'
#' @template ontology
#' @template term_sets
#' @return Character vector of terme IDs.
#' @export
#' @importFrom ontologyIndex get_ancestors
grid_terms <- function(ontology, term_sets) {
	with.ancs <- lapply(term_sets, function(x) get_ancestors(ontology,
		x))
	all.terms <- unique(unlist(with.ancs))
	all.terms[sapply(all.terms, function(term) {
		if (length(intersect(ontology$children[[term]], all.terms)) ==
			0)
			return(TRUE)
		patients.of.each.child <- lapply(intersect(ontology$children[[term]],
			all.terms), function(child) sapply(with.ancs, function(patient.terms) child %in%
			patient.terms))
		patients.of.term <- sapply(with.ancs, function(patient.terms) term %in%
			patient.terms)
		!any(sapply(patients.of.each.child, function(x) identical(x,
			patients.of.term)))
	})]
}

#' Get logical matrix of term annotation for group of cases
#'
#' @template ontology
#' @template term_sets
#' @param all_terms Character vector giving terms to use in annotation.
#' @param remove_unanimous Logical value determining whether to remove terms present in all \code{term_sets}.
#' @param cluster_rows Logical value rows determining whether to use hclust to cluster \code{term_sets}.
#' @param cluster_cols Logical value rows determining whether to use hclust to cluster terms (based on correlation of inclusion in \code{term_sets}). 
#' @return Logical matrix.
#' @export
#' @importFrom stats hclust dist
annotation_grid <- function(ontology, term_sets, all_terms=grid_terms(ontology, term_sets), remove_unanimous=FALSE, cluster_rows=TRUE, cluster_cols=TRUE) {
	with_ancs <- lapply(term_sets, get_ancestors, ontology=ontology)
	mat <- matrix(sapply(all_terms, function(trm) sapply(with_ancs, function(pheno) any(pheno == trm))), ncol=length(all_terms), nrow=length(term_sets), dimnames=list(names(term_sets), ontology$name[all_terms]))
	unan <- apply(mat, 2, all)
	filt_mat <- if (remove_unanimous & (sum(unan) != ncol(mat))) mat[,!unan,drop=FALSE] else mat
	filt_mat[if (cluster_rows & (nrow(filt_mat) > 1)) hclust(dist(filt_mat), method="average")$order else 1:nrow(filt_mat),if (cluster_cols & (ncol(filt_mat) > 1)) hclust(dist(t(filt_mat)), method="average")$order else 1:ncol(filt_mat),drop=FALSE]
}

#' Plot a logical matrix of term annotation
#'
#' @param ... Arguments to be passed to \code{\link{annotation_grid}}.
#' @param on_colour Colour to use to show presence of term.
#' @param off_colour Colour to use to show absence of term.
#' @return Plots heatmap.
#' @export
#' @importFrom paintmap paintmap colour_matrix
plot_annotation_grid <- function(..., on_colour="#FF0000FF", off_colour="#FFFFBFFF") {
	paintmap(colour_matrix(annotation_grid(...), colours=c(off_colour, on_colour)))
}
