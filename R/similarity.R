

term_similarity = function(gl, method = c("kappa", "jaccard", "dice", "overlap"), all = NULL, remove_negative = TRUE) {

    check_pkg("proxyC")
    check_pkg("sparseMatrixStats", bioc = TRUE)

    if(is.null(all)) {
        all = unique(unlist(gl))
    } else {
        gl = lapply(gl, intersect, all)
    }
    gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
    n = length(gl)

    mg = matrix(0, ncol = length(all), nrow = n)
    for(i in seq_len(n)) {
        mg[i, gl[[i]]] = 1
    }
    mg = as(mg, "sparseMatrix")

    method = match.arg(method)[1]
    if(method == "kappa") {
        mat = kappa_dist(mg, remove_negative = remove_negative)
    } else if(method == "overlap") {
        mat = overlap_dist(mg)
    } else {
        mat = proxyC::simil(mg, method = method)
    }

    mat = as.matrix(mat)
    diag(mat) = 1
    rownames(mat) = colnames(mat) = names(gl)
    return(mat)
}


kappa = function(x, y) {
    tab = length(x)
    oab = sum(x == y)/tab
    aab = (sum(x)*sum(y) + sum(!x)*sum(!y))/tab/tab
    k = (oab - aab)/(1 - aab)
    # if(k < 0) k = 0
    return(k)
}

# by rows
kappa_dist = function(m, remove_negative = TRUE) {

    check_pkg("proxyC")
    check_pkg("sparseMatrixStats", bioc = TRUE)

    tab = ncol(m)
    oab = proxyC::simil(m, method = "simple matching")
    m1 = sparseMatrixStats::rowSums2(m)
    m2 = abs(sparseMatrixStats::rowSums2(m - 1))
    aab = (outer(m1, m1) + outer(m2, m2))/tab/tab
    k = (oab - aab)/(1 - aab)
    if(remove_negative) k[k < 0] = 0
    return(k)
}

overlap_dist = function(m) {
    check_pkg("proxyC")
    check_pkg("sparseMatrixStats", bioc = TRUE)
    
    n = sparseMatrixStats::rowSums2(m)
    proxyC::simil(m, method = "dice")*outer(n, n, "+")/2/outer(n, n, pmin)
}

# only for testing
overlap_single = function(x, y) {
    sum(x & y)/min(sum(x), sum(y))
}


#' Jaccard coefficent 
#' 
#' @param lt A list of gene sets.
#' @param weight Weight of genes. The value should be a named vector and it should cover all genes in `lt`.
#' 
#' @export
#' @rdname jaccard
#' @importFrom methods as
#' @examples
#' data(p53_dataset)
#' gs = p53_dataset$gs
#' jaccard(gs[grep("p53", names(gs), ignore.case = TRUE)])
jaccard = function(lt, weight = NULL) {
    n = length(lt)
    nm = names(lt)

    sm = matrix(nrow = n, ncol = n)
    dimnames(sm) = list(nm, nm)
    diag(sm) = 1

    if(n == 1) {
        return(sm)
    }

    universe = unique(unlist(lt, use.names = FALSE))

    if(!is.null(weight)) {
        v = fmatch(universe, names(weight))
        if(any(is.na(v))) {
            stop("`weight` should include all items in `lt`.")
        }
        weight = abs(weight[v])
    }
   
    # convert to integer indicies
    lt = lapply(lt, function(x) {
        v = fmatch(x, universe)
        v[!is.na(v)]
    })
    if(is.null(weight)) {
        ns = sapply(lt, length)
    } else {
        ns = sapply(lt, function(x) sum(weight[x]))
    }

    for(i in 1:(n-1)) {
        for(j in (i+1):n) {
            if(is.null(weight)) {
                v1 = sum(lt[[i]] %fin% lt[[j]])
                v2 = ns[i] + ns[j] - v1
                sm[i, j] = sm[j, i] = v1/v2
            } else {
                ind = lt[[i]][ lt[[i]] %fin% lt[[j]] ]
                v1 = sum(weight[ind])
                v2 = ns[i] + ns[j] - v1
                sm[i, j] = sm[j, i] = v1/v2
            }
        }
    }
    
    return(sm)
}

#' @param min_pct Minimal recurreny of genes in a geneset cluster.
#' @param max_k Maximal number of recurrent genes to use on the plot.
#' @param resolution For controlling the graph clustering method, pass to [`igraph::cluster_louvain()`].
#' @export
#' @rdname jaccard
genesets_similarity_heatmap = function(lt, weight = NULL, min_pct = 0.5, max_k = 5, resolution = 1) {
    sm = jaccard(lt, weight)
    diag(sm) = NA
    
    g = graph_from_adjacency_matrix(sm, weighted = TRUE, mode = "upper", diag = FALSE)
    cl = membership(cluster_louvain(g, resolution = resolution))

    ht = Heatmap(sm, col = c("white", "red"), name = "similarity", 
        border = "black",
        row_split = cl, column_split = cl, show_row_dend = FALSE, 
        show_column_dend = FALSE, row_names_side = "left", 
        row_title = NULL, column_title = NULL)

    ind_list = split(1:nrow(sm), cl)
    ind_list = ind_list[sapply(ind_list, function(ind) {
        p = table(unlist(lt[ind]))/length(ind)
        any(p > min_pct)
    })]

    ind_list = ind_list[sapply(ind_list, length) >= 3]

    if(length(ind_list) == 0) {
        return(draw(ht))
    }

    # recurrent genes
    v_list = lapply(ind_list, function(ind) {
        p = table(unlist(lt[ind]))/length(ind)
        p = sort(p, decreasing = TRUE)
        nm = names(p)[p > min_pct]
        if(length(nm) > max_k) {
            nm = nm[1:max_k]
        }
        df = data.frame(name = nm, p = as.vector(p[nm]))
        if(!is.null(weight)) {
            df$weight = weight[nm]
        }
        df
    })

    if(!is.null(weight)) {
        col_fun_weight = ComplexHeatmap:::default_col(weight, TRUE)
    }

    panel_fun = function(index, nm) {        
        df = v_list[[nm]]
        pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, nrow(df) + 0.5)))
        grid.rect(gp = gpar(fill = "#EEEEEE", col = "#CCCCCC"))
        if(is.null(weight)) {
            grid.rect(0, nrow(df):1, width = df$p, height = 0.6, default.units = "native", just = "left")
        } else {
            grid.rect(0, nrow(df):1, width = df$p, height = 0.6, default.units = "native", just = "left", gp = gpar(fill = col_fun_weight(df$weight)))
        }
        grid.text(df$name, 0, nrow(df):1, gp = gpar(fontsize = 8), hjust = -0.1, default.units = "native")
        grid.annotation_axis(side = "bottom", gp = gpar(fontsize = 6))
        popViewport()
    }

    ht = ht + rowAnnotation(gene = anno_zoom(ind_list, panel_fun = panel_fun, size = unit(sapply(v_list, nrow)*0.5, "cm"), width = unit(3, "cm"),
        link_gp = gpar(fill = "#EEEEEE", col = "#CCCCCC"), internal_line = FALSE, gap = unit(10, "mm")))

    if(is.null(weight)) {
        draw(ht)
    } else {
        draw(ht, heatmap_legend_list = list(Legend(title = "weight", col_fun = col_fun_weight)))
    }
}

