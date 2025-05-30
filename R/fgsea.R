
#' GSEA using fgsea
#' 
#' @param s A numeric vector of gene scores. It should be named with gene IDs.
#' @param gs A list of gene sets. In `fgsea_wrapper()`, genes should have the same ID types as `s`.
#' @param min_size Minimal size of gene sets for analysis.
#' @param max_size Maximal size of gene sets for analysis.
#' @param ... Other argument passed to `fgsea_wrapper()` or further to [`fgsea::fgsea()`].
#' 
#' @rdname fgsea
#' @export
#' @importFrom fgsea fgsea
#' @details
#' Except `fgsea_wrapper()`, gene IDs in `s` in all `fgsea_*()` functions must be EntreZ IDs.
#' 
#' @examples
#' data(p53_dataset)
#' s = p53_dataset$s2n
#' gs = p53_dataset$gs
#' 
#' fgsea_wrapper(s, gs) |> head()
#' 
#' s2 = convert_to_entrez(s)
#' 
#' fgsea_go(s2) |> head()
#' fgsea_kegg(s2) |> head()
#' fgsea_msigdb(s2) |> head()
fgsea_wrapper = function(s, gs, min_size = 5, max_size = 2000, ...) {
    s = sort(s, decreasing = TRUE)
    df = fgsea(gs, s, minSize = min_size, maxSize = max_size, ...)
    df = as.data.frame(df)
    colnames(df) = c("gene_set", "p_value", "p_adjust", "log2_p_err", "es", "nes", "n_gs", "leading_edge")
    df = df[!is.na(df$p_value), , drop = FALSE]
    df = df[order(df$p_value), , drop = FALSE]
    rownames(df) = NULL
    df
}

#' @rdname fgsea
#' @param org_db An `OrgDb` object for the organism. It can be from **org.*.db** packages or downloaded by the **AnnotationHub** package.
#' @param ontology Namespace of GO. Value should be one of "BP", "CC" or "MF".
#' @export
fgsea_go = function(s, org_db = org.Hs.eg.db::org.Hs.eg.db, ontology = "BP", ...) {
    lt = load_go_genesets(org_db, ontology)
   
    df = fgsea_wrapper(s, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname fgsea
#' @param organism See **Details**.
#' @param db A KEGG database. The value can be one of "pathway", "module", "ko", "network", "disease" and "drug".
#' @export
#' @details
#' The value should be set differently for specific `fgsea_*()` functions.
#' 
#' - for `fgsea_kegg()`, the value should be a KEGG organism code, such as "hsa" or "mmu".
#' - for `fgsea_reactome()`, the value should a prefix of the Reactome pathway ID that represents the organism. E.g. "HSA" for human.
#' - for `fgsea_keywords()`, the value can be a organism name, e.g. "human", the latin name or the taxon ID. Please check [`UniProtKeywords::load_keyword_genesets()`].
#' - for `fgsea_phenotype()` and `fgsea_disease()`, the value can only be one of "human", "mouse" and "rat".
#' 
#' All valid values for `fgsea_reactome()` are:
#' 
#' ```
#' c("BTA", "CEL", "CFA", "DRE", "DDI", "DME", "GGA", "HSA", "MMU", 
#'   "MTU", "PFA", "RNO", "SCE", "SPO", "SSC", "XTR")
#' ```
#' 
fgsea_kegg = function(s, organism = "hsa", db = "pathway", ...) {
    lt = load_kegg_genesets(organism, db)

    df = fgsea_wrapper(s, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname fgsea
#' @param collection Collection of the MSigDB gene sets. All possible values can be found via [`list_msigdb_versions()`].
#' @param version Version of the MSigDB database. All possible values can be found via [`list_msigdb_collections()`].
#' @export
fgsea_msigdb = function(s, collection = "h.all", version = "2024.1.Hs", ...) {
    lt = load_msigdb_genesets(collection, version)

    df = fgsea_wrapper(s, lt$gs, ...)
    df
}

#' @rdname fgsea
#' @export
fgsea_reactome = function(s, organism = "HSA", ...) {
    lt = load_reactome_genesets(organism)

    df = fgsea_wrapper(s, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname fgsea
#' @export
fgsea_keywords = function(s, organism = "human", ...) {
    check_pkg("UniProtKeywords", bioc = TRUE)
    gs = UniProtKeywords::load_keyword_genesets(organism)
    
    df = fgsea_wrapper(s, gs, ...)
    df
}

#' @rdname fgsea
#' @export
fgsea_phenotype = function(s, organism = "human", ...) {
    check_pkg("simona", bioc = TRUE)
    lt = load_simona_genesets(simona::ontology_hp, "hp", organism = organism)

    df = fgsea_wrapper(s, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#' @rdname fgsea
#' @export
fgsea_disease = function(s, organism = "human", ...) {
    check_pkg("simona", bioc = TRUE)
    lt = load_simona_genesets(simona::ontology_rdo, "rdo", organism = organism)

    df = fgsea_wrapper(s, lt$gs, ...)
    df$description = lt$names[df$gene_set]
    df
}

#################

#' Native implementation of GSEA
#' 
#' @param expr The complete expression matrix.
#' @param condition The condition labels of samples. The value should be a factor where `level[1]` corresponds to the treatment group and `level[2]` corresponds to the control group.
#' @param gs A list of gene sets. Genes should have the same ID types as `expr`.
#' @param direction Whether the test is one-sided or two-sided.
#' @param power Power added to the absolute value of weight.
#' @param min_size Minimal size of gene sets for analysis.
#' @param max_size Maximal size of gene sets for analysis.
#' @param nperm Number of permutations.
#' @param verbose Whether to print messages?
#' 
#' @rdname native_gsea
#' @export
#' @importFrom matrixStats rowSds
#' @import fastmatch
#' @details
#' The two functions are only for the purpose of studying.
gsea_sample_perm = function(expr, condition, gs, direction = c("std", "pos", "neg"), 
    power = 1, min_size = 5, max_size = 1000, nperm = 1000, verbose = TRUE) {

    if(!is.factor(condition)) {
        stop("`condition` should be a factor.")
    }

    gsea_sample_perm_single = function(expr, condition, gs, power = 1) {
        cmp = levels(condition)
        m1 = expr[, condition == cmp[1]]
        m2 = expr[, condition == cmp[2]]

        s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

        s = sort(s, decreasing = TRUE)
        n = length(s)
        nm = names(s)

        ngs = length(gs)
        es = numeric(ngs)

        for(i in seq_len(ngs)) {
            ind_set = fmatch(gs[[i]], nm)
            l_set = rep(FALSE, n)
            l_set[ind_set] = TRUE

            s_set = abs(s)^power 
            s_set[!l_set] = 0

            l_other = !l_set

            if(direction == "pos") {
                f1 = cumsum(s_set)/sum(s_set)
                f2 = cumsum(l_other)/sum(l_other)
                es[i] = max(f1 - f2)
            } else if(direction == "neg") {
                f1 = cumsum(rev(s_set))/sum(s_set)
                f2 = cumsum(rev(l_other))/sum(l_other)
                es[i] = min(f1 - f2)
            } else {
                f1 = cumsum(s_set)/sum(s_set)
                f2 = cumsum(l_other)/sum(l_other)
                m1 = max(f1 - f2)
                m2 = min(f1 - f2)

                es[i] = max(m1, -m2)*ifelse(m1 > -m2, sign(m1), sign(m2))
            }
        }

        es
    }

    nl = sapply(gs, length)
    gs = gs[nl >= min_size & nl <= max_size]

    es = gsea_sample_perm_single(expr, condition, gs, power)

    ngs = length(gs)
    es_random = matrix(nrow = ngs, ncol = nperm)
    for(i in 1:nperm) {
        es_random[, i] = gsea_sample_perm_single(expr, sample(condition), gs, power)
    }

    if(direction == "std") {
        p = sapply(seq_len(ngs), function(i) {
            if(es[i] > 0) {
                sum(es_random[i, ] > es[i])/sum(es_random[i, ] > 0)
            } else {
                sum(es_random[i, ] < es[i])/sum(es_random[i, ] < 0)
            }
        })
        nes = es/abs(sapply(seq_len(ngs), function(i) {
            if(es[i] > 0) {
                mean(es_random[i, ] > 0)
            } else {
                mean(es_random[i, ] < 0)
            }
        }))
    } else if(direction == "pos") {
        p = sapply(seq_len(ngs), function(i) sum(es_random[i, ] > es[i])/nperm)
        nes = es/rowMeans(es_random)
    } else if(direction == "neg") {
        p = sapply(seq_len(ngs), function(i) sum(es_random[i, ] < es[i])/nperm)
        nes = es/abs(rowMeans(es_random))
    }

    df = data.frame(gene_set = names(gs), p_value = p, p_adjust = p.adjust(p, "BH"), es = es, nes = nes)
    df
}

#' @param s A numeric vector of gene scores.
#' 
#' @rdname native_gsea
#' @export
gsea_gene_perm = function(s, gs, direction = c("std", "pos", "neg"), 
    power = 1, min_size = 5, max_size = 1000, nperm = 1000, verbose = TRUE) {

    s = sort(s, decreasing = TRUE)

    # gene_sets: integer indicies
    direction = match.arg(direction)[1]
    gsea_gene_perm_single = function(s, gs, power = 1, random = FALSE) {
        n = length(s)
        ngs = length(gs)
        es = numeric(ngs)

        for(i in seq_len(ngs)) {
            if(random) {
                ind_set = sample(n, length(gs[[i]]))
            } else {
                ind_set = gs[[i]]
            }

            l_set = logical(n)
            l_set[ ind_set ] = TRUE

            s_set = numeric(n)
            s_set[l_set] = s[l_set]^power
        
            l_other = !l_set

            if(direction == "pos") {
                f1 = cumsum(s_set)/sum(s_set)
                f2 = cumsum(l_other)/sum(l_other)
                es[i] = max(f1 - f2)
            } else if(direction == "neg") {
                f1 = cumsum(rev(s_set))/sum(s_set)
                f2 = cumsum(rev(l_other))/sum(l_other)
                es[i] = min(f1 - f2)
            } else {
                f1 = cumsum(s_set)/sum(s_set)
                f2 = cumsum(l_other)/sum(l_other)
                m1 = max(f1 - f2)
                m2 = min(f1 - f2)

                es[i] = max(m1, -m2)*ifelse(m1 > -m2, sign(m1), sign(m2))
            }
        }

        es
    }

    # absolute weight
    s = abs(s)

    # change genes to integers in the gene sets
    nm = names(s)
    gs = lapply(gs, function(x) {
        ind = fmatch(x, nm)
        ind[!is.na(ind)]
    })
    nl = sapply(gs, length)
    gs = gs[nl >= min_size & nl <= max_size]

    es = gsea_gene_perm_single(s, gs, power)

    # random
    ngs = length(gs)
    es_random = matrix(nrow = ngs, ncol = nperm)
    for(i in 1:nperm) {
        if(verbose) print(i)
        es_random[, i] = gsea_gene_perm_single(s, gs, power, random = TRUE)
    }

    # p-value for different scenarios
    if(direction == "std") {
        p = sapply(seq_len(ngs), function(i) {
            if(es[i] > 0) {
                sum(es_random[i, ] > es[i])/sum(es_random[i, ] > 0)
            } else {
                sum(es_random[i, ] < es[i])/sum(es_random[i, ] < 0)
            }
        })
        nes = es/abs(sapply(seq_len(ngs), function(i) {
            if(es[i] > 0) {
                mean(es_random[i, ] > 0)
            } else {
                mean(es_random[i, ] < 0)
            }
        }))
    } else if(direction == "pos") {
        p = sapply(seq_len(ngs), function(i) sum(es_random[i, ] > es[i])/nperm)
        nes = es/rowMeans(es_random)
    } else if(direction == "neg") {
        p = sapply(seq_len(ngs), function(i) sum(es_random[i, ] < es[i])/nperm)
        nes = es/abs(rowMeans(es_random))
    }

    df = data.frame(gene_set = names(gs), p_value = p, p_adjust = p.adjust(p, "BH"), es = es, nes = nes)
    df
}

