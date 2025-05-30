

#' ID mapping
#' 
#' @param from The "keytype" of the "from" ID type.
#' @param org_db An `OrgDb` object.
#' 
#' @export
#' @importFrom utils getFromNamespace
#' @rdname id_mapping
#' 
#' @details
#' `map_to_entrez()` returns a gene ID mapping from the "from" type to EntreZ IDs, only for protein-coding genes.
map_to_entrez = function (from, org_db = org.Hs.eg.db::org.Hs.eg.db) {

    if(is.character(org_db)) {
        org_db = getFromNamespace(org_db, ns = org_db)
   
    }
    
    x = keys(org_db, keytype = from)
    if("GENETYPE" %in% columns(org_db)) {
        map_tb = AnnotationDbi::select(org_db, keys = x, keytype = from, columns = c("ENTREZID", "GENETYPE"))
        map_tb = map_tb[map_tb$GENETYPE == "protein-coding", , drop = FALSE]
    } else {
        map_tb = AnnotationDbi::select(org_db, keys = x, keytype = from, columns = c("ENTREZID"))
    }
    map_tb = map_tb[!duplicated(map_tb[, 1]), , drop = FALSE]
    
    structure(map_tb[, 2], names = map_tb[, 1])
}


#' @param id A vector of gene IDs
#' @param verbose Whether to print messages?
#' 
#' @export
#' @rdname id_mapping
#' 
#' @details
#' `guess_id_type()` automatically tests the following ID types: ENTREZID, ENSEMBL (genes), REFSEQ and SYMBOL.
guess_id_type = function (id, org_db = org.Hs.eg.db::org.Hs.eg.db, verbose = TRUE) {
    l = grepl("^\\d+$", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENTREZID")
    }
    l = grepl("^ENS.*G", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBL")
    }
    l = grepl("^ENS.*T", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBLTRANS")
    }
    l = grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        return("REFSEQ")
    }

    all_ids = keys(org_db, keytype = "SYMBOL")
    
    l = sample(id, min(100, length(id))) %in% all_ids
    p_match = sum(l)/length(l)
    if (p_match > 0.5) {
        if (verbose) cat("  gene id might be SYMBOL (p = ", sprintf('%.3f', p_match), ")\n")
        return("SYMBOL")
    }
    
    if (verbose) ("  cannot decide which gene id to use.\n")
    return(NULL)
}

guess_id_mapping = function (id, org_db = org.Hs.eg.db::org.Hs.eg.db, verbose = TRUE) {
    col = guess_id_type(id, org_db, verbose = verbose)
    if (is.null(col)) {
        return(NULL)
    }
    if (col == "ENTREZID") {
        return(NULL)
    }
    id_mapping = map_to_entrez(col, org_db)
    l = grepl("^ENS.*(G|T)", id) | grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        fun = local({
            id_mapping = id_mapping
            function(x) {
                x = gsub("\\.\\d+$", "", x)
                id_mapping[x]
            }
        })
        return(fun)
    }
    else {
        return(id_mapping)
    }
}

#' @param x It supports three formats: 1. A vector of gene IDs, 2. A numeric vector with gene IDs as names, 3. A numeric matrix with gene IDs as row names.
#' 
#' @export
#' @rdname id_mapping
#' 
#' @details
#' For `convert_to_entrez()`, if mapping is applied on a character gene ID vector, IDs that cannot be mapped are removed, and for duplicated mappings only 
#' the first one is picked. If mapping is applied on a numeric vector or a matrix, and if there are multiple mappings, the mean is taken as the final value.
convert_to_entrez = function(x, org_db = org.Hs.eg.db::org.Hs.eg.db) {
    if(is.matrix(x)) {
        map = guess_id_mapping(rownames(x), org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        if(is.function(map)) {
            new_rn = map(rownames(x))
        } else {
            new_rn = map[rownames(x)]
        }
        l = is.na(new_rn)

        x = x[!l, , drop = FALSE]
        new_rn = new_rn[!l]

        x2 = do.call(rbind, tapply(1:nrow(x), new_rn, function(ind) {
            colMeans(x[ind, , drop = FALSE])
        }))
        return(x2)

    } else if(is.numeric(x)) {
        map = guess_id_mapping(names(x), org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        x2 = x
        if(is.function(map)) {
            names(x2) = map(names(x))
        } else {
            names(x2) = map[names(x)]
        }
        x2 = x2[!is.na(names(x2))]
        x2 = tapply(x2, names(x2), mean)
        return(structure(as.vector(x2), names = names(x2)))
    } else {
        map = guess_id_mapping(x, org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        if(is.function(map)) {
            x2 = map(x)
        } else {
            x2 = map[x]
        }
        x2 = x2[!is.na(x2)]
        x2 = x2[!duplicated(x2)]
        return(x2)
    }
}

