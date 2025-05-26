
#' Generate gene sets by mapping to orthologues
#' 
#' @param gs A list of gene sets. Genes should be in the EntreZ ID type.
#' @param from Latin name of the "from" organism.
#' @param to Latin name of the "to" organism.
#' 
#' @export
gs_map_to_orthologues = function(gs, from, to) {
    from = gsub(" ", ".", from)
    to = gsub(" ", ".", to)

    check_pkg("Orthology.eg.db", bioc = TRUE)

    keys = keys(Orthology.eg.db::Orthology.eg.db, keytype = from)
    map_df = AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, keys, columns = to, keytype = from)

    map_df[, 1] = as.character(map_df[, 1])
    map_df[, 2] = as.character(map_df[, 2])

    map_vec = structure(map_df[, 2], names = map_df[, 1])

    gs2 = lapply(gs, function(x) {
        x2 = map_vec[x]
        unname(x2[!is.na(x2)])
    })

    gs2[sapply(gs2, length) > 0]
}


CACHE = new.env(parent = emptyenv())

#' Load various gene sets
#' 
#' @param org_db An `OrgDb` object.
#' @param ontology Namespace of the GO. Value should be one of "BP", "CC" and "MF".
#' 
#' @details
#' Genes are all in the EntreZ ID type. All the gene sets are saved in cache.
#' 
#' @export
#' @rdname load_gene_sets
#' @import GO.db
load_go_genesets = function(org_db, ontology) {
    db_info = dbInfo(dbconn(org_db))
    organism = db_info$value[db_info$name == "ORGANISM"]
    if(is.null(CACHE$go_genesets[[organism]][[ontology]])) {
        suppressMessages(tb <- AnnotationDbi::select(org_db, keys = ontology, keytype = "ONTOLOGYALL", columns = c("ENTREZID", "GOALL")))
        gs = split(tb$ENTREZID, tb$GOALL)
        go_gs = lapply(gs, function(x) {
            as.character(unique(x))
        })

        go_names = Term(GOTERM)[names(go_gs)]

        CACHE$go_genesets[[organism]][[ontology]]$go_names = go_names
        CACHE$go_genesets[[organism]][[ontology]]$go_gs = go_gs
    } else {
        go_names = CACHE$go_genesets[[organism]][[ontology]]$go_names
        go_gs = CACHE$go_genesets[[organism]][[ontology]]$go_gs
    }

    list(names = go_names, gs = go_gs)
}

#' @param organism See **Details**.
#' @param db A KEGG database. The value can be one of "pathway", "module", "ko", "network", "disease" and "drug".
#' @export
#' @details
#' The value should be set differently for specific `fgsea_*()` functions.
#' 
#' - for `load_kegg_genesets()`, the value should be a KEGG organism code, such as "hsa" or "mmu".
#' - for `load_reactome_genesets()`, the value should a prefix of the Reactome pathway ID that represents the organism. E.g. "HSA" for human.
#' - for `load_simona_genesets()`, the value can only be one of "human", "mouse" and "rat".
#' 
#' @rdname load_gene_sets
load_kegg_genesets = function(organism, db) {
    if(is.null(CACHE$kegg_pathways[[db]][[organism]])) {
        if(db == "pathway") {
            gs_names = read.table(url(paste0("https://rest.kegg.jp/list/", db, "/", organism)), sep = "\t", quote = "")
        } else {
            gs_names = read.table(url(paste0("https://rest.kegg.jp/list/", db)), sep = "\t", quote = "")
        }
        gs_names = structure(gs_names[, 2], names = gs_names[, 1])

        gs_genes = read.table(url(paste0("https://rest.kegg.jp/link/", db, "/", organism)), sep = "\t", quote = "")
        gs_genes[, 1] = gsub("^\\w+:", "", gs_genes[, 1])
        gs_genes[, 2] = gsub("^\\w+:", "", gs_genes[, 2])
        if(db == "module") {
            gs_genes[, 2] = gsub("^[a-zA-Z]+_", "", gs_genes[, 2])
        }

        gs = split(gs_genes[, 1], gs_genes[, 2])
        gs = lapply(gs, function(x) {
            as.character(unique(x))
        })

        CACHE$kegg_pathways[[db]][[organism]]$gs_names = gs_names
        CACHE$kegg_pathways[[db]][[organism]]$gs = gs
    } else {
        gs_names = CACHE$kegg_pathways[[db]][[organism]]$gs_names
        gs = CACHE$kegg_pathways[[db]][[organism]]$gs
    }

    list(names = gs_names, gs = gs)
}

#' @param collection Collection of the MSigDB gene sets. All possible values can be found via [`list_msigdb_versions()`].
#' @param version Version of the MSigDB database. All possible values can be found via [`list_msigdb_collections()`].
#' @export
#' @rdname load_gene_sets
load_msigdb_genesets = function(collection, version) {
    dataset = paste0(version, "/", collection)
    if(is.null(CACHE$msigdb_genesets[[dataset]])) {
        gs = get_msigdb(version = version, collection = collection)
        CACHE$msigdb_genesets[[dataset]] = gs
    } else {
        gs = CACHE$msigdb_genesets[[dataset]]
    }

    list(gs = gs)
}

#' @export
#' @rdname load_gene_sets
load_reactome_genesets = function(organism) {
    if(is.null(CACHE$reactome_pathways)) {

        check_pkg("reactome.db", bioc = TRUE)

        pathway_names = unlist(as.list(reactome.db::reactomePATHID2NAME))
        pathway_gs = as.list(reactome.db::reactomePATHID2EXTID)

        cn = intersect(names(pathway_names), names(pathway_gs))
        pathway_names = pathway_names[cn]
        pathway_gs = pathway_gs[cn]

        CACHE$reactome_pathways$pathway_gs = pathway_gs
        CACHE$reactome_pathways$pathway_names = pathway_names
    } else {
        pathway_gs = CACHE$reactome_pathways$pathway_gs
        pathway_names = CACHE$reactome_pathways$pathway_names
    }

    l = grepl(paste0("^R-", organism, "-"), names(pathway_names))
    list(names = pathway_names[l], gs = pathway_gs[l])
}

#' @param fun Use `ontology_hp` and `ontology_rdo`.
#' @param ontology Just a tag, please use "hp" and "rdo".
#' 
#' @details
#' The human phenotype gene sets and disease gene sets are supported by the **simona** package.
#' 
#' @export
#' @rdname load_gene_sets
#' @importFrom S4Vectors mcols
load_simona_genesets = function(fun, ontology, organism = "human") {
    if(is.null(CACHE$simona_genesets[[ontology]][[organism]])) {

        check_pkg("simona", bioc = TRUE)

        dag = fun(organism)
        mc = mcols(dag)
        gs_names = structure(mc$name, names = rownames(mc))
        all_terms = simona::dag_all_terms(dag)
        gs = simona::term_annotations(dag, all_terms)
        gs_names = gs_names[all_terms]

        org_db = switch(organism,
            human = org.Hs.eg.db::org.Hs.eg.db,
            mouse = {check_pkg("org.Ms.eg.db", bioc = TRUE); org.Ms.eg.db::org.Ms.eg.db},
            rat   = {check_pkg("org.Rn.eg.db", bioc = TRUE); org.Rn.eg.db::org.Rn.eg.db},
            stop("not supported")
        )

        x = keys(org_db, keytype = "SYMBOL")
        suppressMessages(map_tb <- AnnotationDbi::select(org_db, keys = x, keytype = "SYMBOL", columns = c("ENTREZID", "GENETYPE")))
        map_tb = map_tb[map_tb$GENETYPE == "protein-coding", , drop = FALSE]
        map_tb = map_tb[!duplicated(map_tb[, 1]), , drop = FALSE]
        map = structure(map_tb[, 2], names = map_tb[, 1])

        gs = lapply(gs, function(x) {
            x2 = map[x]
            unique(x2[!is.na(x2)])
        })

        n = sapply(gs, length)
        gs = gs[n > 0]
        gs_names = gs_names[n > 0]

        CACHE$simona_genesets[[ontology]][[organism]]$gs = gs
        CACHE$simona_genesets[[ontology]][[organism]]$gs_names = gs_names
    } else {
        gs = CACHE$simona_genesets[[ontology]][[organism]]$gs
        gs_names = CACHE$simona_genesets[[ontology]][[organism]]$gs_names
    }

    list(gs = gs, names = gs_names)
}
