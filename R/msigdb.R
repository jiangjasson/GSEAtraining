
msigdb_base_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
msigdb_env = new.env(parent = emptyenv())
msigdb_env$all_versions = NULL
msigdb_env$file_list = list()
msigdb_env$all_collections = list()
msigdb_env$gene_sets = list()


#' Retrieve gene sets from MSigDB
#' 
#' @rdname msigdb
#' @export
#' @examples
#' list_msigdb_versions()
#' list_msigdb_collections("2024.1.Hs")
list_msigdb_versions = function() {
    if(is.null(msigdb_env$all_versions)) {
        all_versions = get_file_list(msigdb_base_url)
        msigdb_env$all_versions = all_versions
    } else {
        all_versions = msigdb_env$all_versions
    }
    all_versions
}

#' @importFrom utils menu
choose_msigdb_version = function() {
    all_versions = list_msigdb_versions()
    ind = menu(all_versions, title = "Choose an MSigDB version:")
    all_versions[ind]
}

#' @rdname msigdb
#' 
#' @param version The MSigDB version. The value can be obtained by `list_msigdb_versions()`. 
#'         If this argument is not specified, a menu will be printed to let users to select.
#' @export
list_msigdb_collections = function(version = NULL) {
    if(is.null(version)) {
        version = choose_msigdb_version()
    }

    if(is.null(msigdb_env$file_list[[version]])) {
        files = get_file_list(paste0(msigdb_base_url, "/", version))
        files = files[grep("\\.gmt$", files)]
        files = files[!grepl("^msigdb", files)]
        msigdb_env$file_list[[version]] = files
    } else {
        files = msigdb_env$file_list[[version]]
    }
    collections = unique(gsub(paste0(".v", version, ".*$"), "", files))
    msigdb_env$all_collections[[version]] = collections
    collections
}

#' @rdname msigdb
#' @export
choose_msigdb_collection = function(version) {
    
    collections = list_msigdb_collections(version)
    ind = menu(collections, title = paste0("Choose a gene set collection for version ", version, ":"))
    collections[ind]
}

get_file_list = function(url) {
    con = url(url)
    on.exit(close(con))
    ln = readLines(con)
    ind = grep("^<tr><td", ln)
    rows = ln[ind]
    rows = gsub("</td><td[^>]*>", ";", rows)
    rows = gsub("<.*?>", "", rows)
    files = sapply(strsplit(rows, ";"), function(x) x[2])[-1]
    gsub("/", "", files)
}


#' @rdname msigdb
#' 
#' @param collection The gene set collection. The values can be obtained by `list_msigdb_collections()`.
#' @param gene_id_type One of "entrez" and "symbols".
#' @param as_table Whether to return a list or a table.
#' @export
get_msigdb = function(version = choose_msigdb_version(), 
    collection = choose_msigdb_collection(version),
    gene_id_type = c("entrez", "symbols"), as_table = FALSE) {
    
    version = force(version)
    gene_id_type = match.arg(gene_id_type)

    if(is.null(msigdb_env$all_versions)) {
        list_msigdb_versions()
    }
    
    if(!version %in% msigdb_env$all_versions) {
        i = grep(version, msigdb_env$all_versions, ignore.case = TRUE)
        if(length(i) == 1) {
            version = msigdb_env$all_versions[i]
        } else {
            message(paste0("Cannot find version '", version, "', please select a valid value."))
            version = choose_msigdb_version()
        }
    }

    collection = force(collection)
    if(is.null(msigdb_env$all_collections[[version]])) {
        list_msigdb_collections(version)
    }
    if(!collection %in% msigdb_env$all_collections[[version]]) {
        i = grep(collection, msigdb_env$all_collections[[version]], ignore.case = TRUE)
        if(length(i) == 1) {
            collection = msigdb_env$all_collections[[version]][i]
        } else {
            message(paste0("Cannot find collection '", collection, "', please select a valid value."))
            collection = choose_msigdb_collection(version)
        }
    }
    
    url = paste0(msigdb_base_url, "/", version, "/", collection, ".v", version, ".", gene_id_type, ".gmt")
    basename = basename(url)
    if(is.null(msigdb_env$gene_sets[[basename]])) {
        con = url(url)
        on.exit(close(con))
        ln = readLines(con)
        ln = strsplit(ln, "\t")
        gs = lapply(ln, function(x) x[-(1:2)])
        names(gs) = sapply(ln, function(x) x[1])
        msigdb_env$gene_sets[[basename]] = gs
    } else {
        gs = msigdb_env$gene_sets[[basename]]
    }
   
    if(as_table) {
        df = data.frame(gene_set = rep(names(gs), times = sapply(gs, length)),
                        gene = unlist(gs))
        rownames(df) = NULL
        df
    } else {
        gs
    }
}
