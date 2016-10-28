#' Get an adjacency matrix for a set of ontological terms
#' 
#' @template ontology
#' @template terms
#' @return A logical matrix representing the adjacency matrix of \code{terms} based on the directed acyclic graph of \code{ontology}. A \code{TRUE} entry means the term correspnding to the column is a parent of the row term within \code{terms}.
#' @seealso \code{\link{get_adjacency_matrix}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_pseudo_adjacency_matrix(hpo, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
#' @importFrom ontologyIndex minimal_set
#' @importFrom stats setNames
get_pseudo_adjacency_matrix <- function(ontology, terms) structure(
	t(
		sapply(
			setNames(terms, terms),
			function(term) "%in%"(
				setNames(terms, terms),
				minimal_set(	
					ontology, 
					setdiff(
						intersect(terms, ontology$ancestors[[term]]),
						term
					)
				)
			)
		)
	),
	dimnames=rep(
		list(terms),
		2
	)
)

remove_uninformative_once <- function(ontology, with.ancs, all.terms) {
	mat <- get_pseudo_adjacency_matrix(ontology, all.terms)
	has_term <- sapply(with.ancs, function(x) all.terms %in% x)

	all.terms[sapply(seq(nrow(mat)), function(term_index) {
		if (sum(mat[,term_index]) == 0) TRUE
		else !all(sapply(split(has_term[mat[,term_index],,drop=FALSE], seq(sum(mat[,term_index]))), identical, unname(has_term[term_index,])))
	})]
}

#' Remove uninformative terms from union of all terms in set of annotations
#' 
#' For a set of ontological annotation sets, remove terms annotated to the same objects as all their children. Useful for selecting terms for summarising a set of annotation sets, as it can lead to a significant reduction in the number of terms.
#'
#' @template ontology
#' @template term_sets
#' @return Character vector of terms
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' remove_uninformative_terms(hpo, list(Patient1=c("HP:0001873","HP:0000118")))
#' @export
#' @importFrom ontologyIndex get_ancestors
remove_uninformative_terms <- function(ontology, term_sets) {
	with.ancs <- lapply(term_sets, get_ancestors, ontology=ontology)
	all.terms <- unique(unlist(use.names=FALSE, with.ancs))
	repeat {
		new_terms <- remove_uninformative_once(ontology, with.ancs, all.terms)
		if (setequal(new_terms, all.terms))
			break
		else
			all.terms <- new_terms
	}
	new_terms
}

#' Get an adjacency matrix for a set of ontological terms
#' 
#' @template ontology
#' @template terms
#' @return A logical matrix representing the adjacency matrix of \code{terms} based on the directed acyclic graph of \code{ontology}. A \code{TRUE} entry means the term correspnding to the column is a parent of the row term in \code{ontology}.
#' @seealso \code{\link{get_pseudo_adjacency_matrix}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_adjacency_matrix(hpo, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
get_adjacency_matrix <- function(ontology, terms) {
	names(terms) <- terms
	adj.mat <- sapply(
		terms,
		function(term) terms %in% ontology$parents[[term]]
	)
	rownames(adj.mat) <- colnames(adj.mat)
	t(adj.mat)
}

