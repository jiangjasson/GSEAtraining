
#' @rdname topology_kegg
#' @export
#' @importFrom utils download.file
retrieve_all_kegg_graphs = function(organism = "hsa") {

	check_pkg("KEGGgraph", bioc = TRUE)

	tb = read.table(url(paste0("https://rest.kegg.jp/list/pathway/", organism)), sep = "\t")
	all_pathway_ids = tb[, 1]

	pl = list()
	tmp = tempfile()
	for(i in seq_along(all_pathway_ids)) {
		cat("download pathway kgml for", all_pathway_ids[i], "\n")
		url = paste0("https://rest.kegg.jp/get/", all_pathway_ids[i], "/kgml")
		download.file(url, destfile = tmp)
		pl[[ all_pathway_ids[i] ]] = KEGGgraph::parseKGML2DataFrame(tmp)
	}

	pl2 = lapply(pl, function(x) {

		if(nrow(x) == 0) {
			return(NULL)
		}
		x = x[, 1:2]
		x[[1]] = gsub("^\\w+:", "", x[[1]])
		x[[2]] = gsub("^\\w+:", "", x[[2]])
		graph_from_data_frame(x)
	})
	pl2 = pl2[!sapply(pl2, is.null)]

	attr(pl2, "timestamp") = Sys.Date()
	pl2
}


# pl = retrieve_all_kegg_graphs("hsa")
# saveRDS(pl, file = "~/project/GSEAtraining/inst/extdata/kegg_hsa_pathway_graphs.rds", compress = "xz")
# pl = retrieve_all_kegg_graphs("mmu")
# saveRDS(pl, file = "~/project/GSEAtraining/inst/extdata/kegg_mmu_pathway_graphs.rds", compress = "xz")
# pl = retrieve_all_kegg_graphs("rno")
# saveRDS(pl, file = "~/project/GSEAtraining/inst/extdata/kegg_rno_pathway_graphs.rds", compress = "xz")

#' Centrality-based KEGG enrichment analysis
#' 
#' @param genes A vector of genes.
#' @param universe A vector of universe genes. If it is not specified, the total genes in the KEGG pathways are used.
#' @param centrality Centrality method. The value should be a function which accepts an `igraph` object and returns a vector of centrality values.
#' @param organism KEGG organism code.
#' @param pl A list of KEGG pathways as `igraph` objects. If the organism is one of "hsa", "mmu" or "rno", the corresponding `pl` is already
#'      generated and will be loaded automatically. For other organisms, use `retrieve_all_kegg_graphs()` to generate one.
#' @param perm Number of permutations.
#' @param min_hits Minimal number of the overlapping genes in `genes` and pathways.
#' @param min_size Minimal number of genes in pathways.
#' @param max_size Maximal number of genes in pathways.
#' @param verbose Whether to print messages?
#' 
#' @details
#' The following lists several useful centrality measures, written as functions. They can be assigned to the `centrality` argument.
#' 
#' - in-degree: `function(g) igraph::degree(g, mode = "in")`
#' - out-degree: `function(g) igraph::degree(g, mode = "out")`
#' - betweenness: `igraph::betweenness`
#' - page rank: `igraph::page_rank`
#' - in-reach: `function(g) CePa::reach(g, mode = "in")`
#' - out-reach: `function(g) CePa::reach(g, mode = "out")`
#' 
#' And many more in the **igraph** and **CePa** packages.
#' 
#' `ora_kegg_topology()` is the ORA-extension.
#' 
#' @export
#' @import igraph
#' @rdname topology_kegg
ora_kegg_topology = function(genes, universe = NULL, centrality = igraph::degree, organism = "hsa", pl = NULL,
    perm = 1000, min_hits = 3, min_size = 5, max_size = 2500, verbose = TRUE) {

    if(organism %in% c("hsa", "mmu", "rno")) {
        pl = readRDS(system.file("extdata", paste0("kegg_", organism, "_pathway_graphs.rds"), package = "GSEAtutorial"))
    } else {
        if(is.null(pl)) {
            stop("please use `retrieve_all_kegg_graphs()` to generate the pathway data and assign it to `pl` argument.")
        }
    }

    if(is.null(universe)) {
    	universe = unique(unlist(lapply(pl, function(x) V(x)$name), use.names = FALSE))
    }

    cen = lapply(pl, function(g) {
        centrality(g)
    })

    cen = lapply(cen, function(x) {
        cn = intersect(names(x), universe)
        x[cn]
    })

    nl = sapply(cen, length)
    cen = cen[nl >= min_size & nl <= max_size]
    nl = nl[nl >= min_size & nl <= max_size]
    
    stat = sapply(cen, function(x) {
        cn = intersect(genes, names(x))
        c(k = length(cn), s = sum(x[cn]))
    })
    k = stat["k", ]
    s = stat["s", ]

    l = k >= min_hits
    k = k[l]
    s = s[l]
    cen = cen[l]
    nl = nl[l]

    ngs = length(cen)

    stat_random = matrix(nrow = ngs, ncol = perm)
    for(i in 1:perm) {
        if(verbose) {
            cat(strrep("\b", 100)); cat(strrep(" ", 100)); cat(strrep("\b", 100)); cat(i, "/", perm)
        }
        genes_random = sample(universe, length(genes))
        stat_random[, i] = sapply(cen, function(x) {
            cn = intersect(genes_random, names(x))
            sum(x[cn])
        })
    }
    if(verbose) cat("\n")

    p = sapply(1:ngs, function(i) {
        l = stat_random[i, ] > s[i]
        l[is.na(l)] = FALSE
        sum(l)/perm
    })

    log2fe = log2(s/rowMeans(stat_random, na.rm = TRUE))
    z = (s - rowMeans(stat_random, na.rm = TRUE))/rowSds(stat_random, na.rm = TRUE)

    df = data.frame(pathway = names(cen), k = k, s = s, n_gs = nl, log2fe = log2fe, z_score = z, p_value = p, p_adjust = p.adjust(p, "BH"))

    tb = read.table(url(paste0("https://rest.kegg.jp/list/pathway/", organism)), sep = "\t")
    pathway_names = structure(tb[, 2], names = tb[, 1])
    df$description = pathway_names[df$pathway]

    df = df[order(df$p_value), ]
    rownames(df) = NULL
    df
}

#' @rdname topology_kegg
#' @export
#' @param s A numeric vector of gene scores.
#' 
#' @details
#' `gsea_kegg_topology()` is the GSEA extension. The geneset-level score is calculated as `mean(abs(s*w))`.
#' 
#' @importFrom utils read.table
gsea_kegg_topology = function(s, centrality = igraph::degree, organism = "hsa", pl = NULL,
    min_size = 5, max_size = 2500, verbose = TRUE) {

    if(organism %in% c("hsa", "mmu", "rno")) {
        pl = readRDS(system.file("extdata", paste0("kegg_", organism, "_pathway_graphs.rds"), package = "GSEAtutorial"))
    } else {
        if(is.null(pl)) {
            stop("please use `retrieve_all_kegg_graphs()` to generate the pathway data and assign it to `pl` argument.")
        }
    }

    cen = lapply(pl, function(g) {
        centrality(g)
    })

    universe = names(s)
    cen = lapply(cen, function(x) {
        cn = intersect(names(x), universe)
        x[cn]
    })

    nl = sapply(cen, length)
    cen = cen[nl >= min_size & nl <= max_size]
    nl = nl[nl >= min_size & nl <= max_size]
    
    stat = sapply(cen, function(x) {
        mean(abs(s[names(x)]) * x)
    })
    
    ngs = length(cen)

    stat_random = matrix(nrow = ngs, ncol = 1000)
    for(i in 1:1000) {
        if(verbose) {
            cat(strrep("\b", 100)); cat(strrep(" ", 100)); cat(strrep("\b", 100)); cat(i, "/", 1000)
        }
        stat_random[, i] = sapply(cen, function(x) {
            s_random = sample(s, length(x))
            mean(abs(s_random)*x)
        })
    }
    if(verbose) cat("\n")

    p = sapply(1:ngs, function(i) {
        l = stat_random[i, ] > stat[i]
        l[is.na(l)] = FALSE
        sum(l)/1000
    })

    log2fe = log2(stat/rowMeans(stat_random, na.rm = TRUE))
    z = (stat - rowMeans(stat_random, na.rm = TRUE))/rowSds(stat_random, na.rm = TRUE)

    df = data.frame(pathway = names(cen), s = stat, n_gs = nl, log2fe = log2fe, z_score = z, p_value = p, p_adjust = p.adjust(p, "BH"))

    tb = read.table(url(paste0("https://rest.kegg.jp/list/pathway/", organism)), sep = "\t")
    pathway_names = structure(tb[, 2], names = tb[, 1])
    df$description = pathway_names[df$pathway]

    df = df[order(df$p_value), ]
    rownames(df) = NULL
    df
}
