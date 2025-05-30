
#' Get genes and GO terms from OrgDb oject
#' 
#' @param org_db An `OrgDb` object.
#' @param ontology "BP", "CC" or "MF".
#' @param gene_id_type Which gene ID type to use?
#' @param as_table Whether to return a list or a two-column data frame?
#' 
#' @export
#' @rdname orgdb
#' @examples
#' require(org.Hs.eg.db)
#' gs = get_GO_gene_sets_from_orgdb(org.Hs.eg.db)
#' gs[1:2]
get_GO_gene_sets_from_orgdb = function(org_db, ontology = "BP", 
    gene_id_type = c("ENTREZID", "SYMBOL", "ENSEMBL"), as_table = FALSE) {
    
    ontology = toupper(ontology)

    gene_id_type = match.arg(gene_id_type)[1]

    tb = AnnotationDbi::select(org_db, keys = ontology, keytype = "ONTOLOGYALL", columns = c(gene_id_type, "GOALL"))
    tb = tb[, 2:3]
    tb = unique(tb)
    if(as_table) {
        tb
    } else {
        split(tb[[1]], tb[[2]])
    }
}


#' @param n Number of random genes.
#' 
#' @export
#' @import AnnotationDbi
#' @rdname orgdb
#' @examples
#' random_genes(org.Hs.eg.db) |> head()
random_genes = function(org_db, n = 1000, gene_id_type = c("ENTREZID", "SYMBOL", "ENSEMBL")) {
    gene_id_type = match.arg(gene_id_type)[1]
    if("GENETYPE" %in% columns(org_db)) {
        gi = AnnotationDbi::select(org_db, key = "protein-coding", keytype = "GENETYPE", column = gene_id_type)[, 2]
    } else {
        gi = keys(org_db, keytype = gene_id_type)
    }
    sample(gi, min(n, length(gi)))
}

#' @export
#' @rdname orgdb
#' @examples
#' protein_coding_genes(org.Hs.eg.db) |> head()
protein_coding_genes = function(org_db, gene_id_type = c("ENTREZID", "SYMBOL", "ENSEMBL")) {
    gene_id_type = match.arg(gene_id_type)[1]
    AnnotationDbi::select(org_db, key = "protein-coding", keytype = "GENETYPE", column = gene_id_type)[, 2]
}
