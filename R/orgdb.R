
#' Get GO gene sets from OrgDb object
#' 
#' @param org_db An `OrgDb` object.
#' @param ontology "BP", "CC" or "MF".
#' @param gene_id_type Which gene ID type to use.
#' @param as_table Whether to return a list or a two-column data frame?
#' 
#' @export
get_GO_gene_sets_from_orgdb = function(org_db, ontology = "BP", 
    gene_id_type = c("ENTREZID", "SYMBOL", "ENSEMBL"), as_table = FALSE) {
    
    ontology = toupper(ontology)

    gene_id_type = match.arg(gene_id_type)[1]

    tb = AnnotationDbi::select(org_db, keys = keys(org_db, gene_id_type), columns = c("GOALL", "ONTOLOGYALL"), keytype = gene_id_type)
    tb = tb[tb$ONTOLOGYALL %in% ontology, , drop = FALSE]
    tb = tb[!is.na(tb$GOALL), c(gene_id_type, "GOALL"), drop = FALSE]
    tb = unique(tb)
    if(as_table) {
        tb
    } else {
        split(tb[, gene_id_type], tb$GOALL)
    }
}

#' Random genes
#' 
#' @param org_db An `OrgDb` object.
#' @param n Number of random genes.
#' @param keytype Keytype of the genes in the `OrgDb` database.
#' 
#' @export
#' @import AnnotationDbi
random_genes = function(org_db, n = 1000, keytype = "ENTREZID") {
    if("GENETYPE" %in% columns(org_db)) {
        gi = AnnotationDbi::select(org_db, key = "protein-coding", keytype = "GENETYPE", column = keytype)[, 2]
    } else {
        gi = keys(org_db, keytype = keytype)
    }
    sample(gi, min(n, length(gi)))
}

protein_coding_genes = function(org_db, keytype = "ENTREZID") {
    AnnotationDbi::select(org_db, key = "protein-coding", keytype = "GENETYPE", column = keytype)[, 2]
}
