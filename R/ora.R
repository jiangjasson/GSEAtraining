

#' Over-representation analysis
#' 
#' @param genes A vector of genes.
#' @param gs A list of gene sets. Genes should have the smae ID type as in `genes`.
#' @param universe A vector of background genes.
#' @param min_hits Minimal number of overlapping genes in `genes` and gene sets.
#' @param min_size Minimal number of genes in gene sets.
#' @param max_size Maximal number of genes in gene sets.
#' 
#' @importFrom stats phyper
#' @importFrom stats p.adjust
#' @export
#' @import fastmatch
#' @rdname ora
ora = function(genes, gs, universe = NULL, min_hits = 3, min_size = 5, max_size = 2500) {

    if(is.null(universe)) {
        universe = unique(unlist(gene_sets, use.names = FALSE))
    } else {
        universe = unique(universe)
    }

    nl = sapply(gene_sets, length)
    gene_sets = gene_sets[nl >= min_size & nl <= max_size]
    
    gs_name = names(gene_sets)
    
    genes = fmatch(genes, universe); genes = unique(genes[!is.na(genes)])
    gene_sets = lapply(gene_sets, function(x) {
        v = fmatch(x, universe)
        unique(v[!is.na(v)])
    })

    n_universe = length(universe)
    k = length(genes)
    
    x = sapply(gene_sets, function(x) sum(x %fin% genes))
    m = sapply(gene_sets, length)
    if(all(m == 0)) {
        stop("Gene IDs in `genes` and `gene_sets` do not match.")
    }
    n = n_universe - m
    
    p = phyper(x - 1, m, n, k, lower.tail = FALSE)

    df = data.frame(
        gene_set = gs_name,
        n_hits = x,
        n_genes = k,
        n_gs = m,
        n_total = n_universe,
        log2fe = log2(x*n_universe/k/m),
        z_score = (x - k*m/n_universe)/sqrt(m*k/n_universe*(n_universe-k)/n_universe*(n_universe-m)/(n_universe-1)),
        p_value = p
    )
    df$p_adjust = p.adjust(df$p_value, "BH")
    rownames(df) = NULL
    df[order(df$p_adjust, df$p_value), ,drop = FALSE]
}


#' @rdname ora
#' @param ... Passed to `ora()`.
#' @export
#' @param org_db An `OrgDb` object for the organism. It can be from **org.*.db** packages or downloaded by the **AnnotationHub** package.
#' @param ontology Namespace of GO. Value should be one of "BP", "CC" or "MF".
ora_go = function(genes, org_db = org.Hs.eg.db::org.Hs.eg.db, ontology = "BP", ...) {
    lt = load_go_genesets(org_db, ontology)
   
    df = ora(genes, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname ora
#' @param organism See **Details**.
#' @param db A KEGG database. The value can be one of "pathway", "module", "ko", "network", "disease" and "drug".
#' @export
#' @details
#' The value should be set differently for specific `ora_*()` functions.
#' 
#' - for `ora_kegg()`, the value should be a KEGG organism code, such as "hsa" or "mmu".
#' - for `ora_reactome()`, the value should a prefix of the Reactome pathway ID that represents the organism. E.g. "HSA" for human.
#' - for `ora_keywords()`, the value can be a organism name, e.g. "human", the latin name or the taxon ID.
#' - for `ora_phenotype()` and `fgsea_disease()`, the value can only be one of "human", "mouse" and "rat".
#' 
ora_kegg = function(genes, organism = "hsa", db = "pathway", ...) {
    lt = load_kegg_genesets(organism, db)

    df = ora(genes, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname ora
#' @param collection Collection of the MSigDB gene sets. All possible values can be found via [`list_msigdb_versions()`].
#' @param version Version of the MSigDB database. All possible values can be found via [`list_msigdb_collections()`].
#' @export
ora_msigdb = function(genes, collection = "h.all", version = "2024.1.Hs", ...) {
    lt = load_msigdb_genesets(collection, version)

    df = ora(genes, lt$gs, ...)
    df
}

#' @rdname ora
#' @export
ora_reactome = function(genes, organism = "HSA", ...) {
    lt = load_reactome_genesets(organism)

    df = ora(genes, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname ora
#' @export
ora_keywords = function(genes, organism = "human", ...) {
    check_pkg("UniProtKeywords", bioc = TRUE)
    gs = UniProtKeywords::load_keyword_genesets(organism)
    df = ora(genes, gs, ...)
    df
}


#' @rdname ora
#' @export
ora_phenotype = function(genes, organism = "human", ...) {
    check_pkg("simona", bioc = TRUE)
    lt = load_simona_genesets(simona::ontology_hp, "hp", organism = organism)

    df = ora(genes, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname ora
#' @export
ora_disease = function(genes, organism = "human", ...) {
    check_pkg("simona", bioc = TRUE)
    lt = load_simona_genesets(simona::ontology_rdo, "rdo", organism = organism)

    df = ora(genes, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}
