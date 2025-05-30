

#' Single-sample GSEA transformation
#' 
#' @param mat A numeric matrix.
#' @param gs A list of gene sets. Genes should have the same ID type as in `mat`.
#' @param scale Whether to apply z-score scaling on rows?
#' @param set_level A function to calculate the geneset-level score. This function should accept two arguments: a vector of scores for all genes, and a corresponding logical vector repressenting which genes are in the current gene set.
#' 
#' @export
#' @examples
#' data(p53_dataset)
#' expr = p53_dataset$expr
#' gs = p53_dataset$gs
#' mat_gs = ss_gsea(expr, gs)
#' head(mat_gs)
ss_gsea = function(mat, gs, scale = TRUE, set_level = set_level_mean) {
	if(scale) {
		mat = t(scale(t(mat)))
	}

	mat_gs = matrix(nrow = length(gs), ncol = ncol(mat))
	rownames(mat_gs) = names(gs)
	colnames(mat_gs) = colnames(mat)
	for(i in 1:ncol(mat)) {
		for(j in seq_along(gs)) {
	   		mat_gs[j, i] = set_level(mat[, i], rownames(mat) %in% gs[[j]])
		}
	}

	mat_gs
}
