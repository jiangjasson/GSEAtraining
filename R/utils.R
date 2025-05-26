

#' Convert gene sets between list and a data frame
#' 
#' @param list A list of genes.
#' @rdname list_to_dataframe
#' @export
list_to_dataframe = function(list) {
    df = data.frame(gene_set = rep(names(list), times = sapply(list, length)),
               gene = unlist(list))
    rownames(df) = NULL
    df
}


#' @param df A two-column data frame
#' 
#' @details
#' In `gs_dataframe_to_list()`, which column contains genes and which column contains gene sets
#' are automatically checked by the number of genes and gene sets. Basically number of genes
#' should be larger than the number of gene sets.
#' 
#' @rdname list_to_dataframe
#' @export
dataframe_to_list = function(df) {
    n1 = length(unique(df[, 1]))
    n2 = length(unique(df[, 2]))
    if(n1 < n2) {
        split(df[, 2], df[, 1])
    } else {
        split(df[, 1], df[, 2])
    }
}



#' Wrap long text into several lines
#' 
#' @param x A vector of sentences.
#' @param width Maximal number of chararacters per line.
#' @export
wrap_text = function(x, width = 60) {
    sapply(x, function(txt) {
        paste(strwrap(txt, width = width), collapse = "\n")
    })
}

#' @importFrom utils install.packages
check_pkg = function(pkg, bioc = FALSE) {
    if(requireNamespace(pkg, quietly = TRUE)) {
        return(NULL)
    } else {

        if(!interactive()) {
            if(bioc) {
                stop(paste0("You need to manually install package '", pkg, "' from Bioconductor."))
            } else {
                stop(paste0("You need to manually install package '", pkg, "' from CRAN."))
            }
        }

        if(bioc) {
            answer = readline(paste0("Package '", pkg, "' is required but not installed. Do you want to install it from Bioconductor? [y|n] "))
        } else {
            answer = readline(paste0("Package '", pkg, "' is required but not installed. Do you want to install it from CRAN? [y|n] "))
        }

        if(bioc) {
            if(tolower(answer) %in% c("y", "yes")) {
                if(!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager")
                }
                BiocManager::install(pkg)
            } else {
                stop(paste0("You need to manually install package '", pkg, "' from Bioconductor."))
            }
        } else {
            if(tolower(answer) %in% c("y", "yes")) {
                install.packages(pkg)
            } else {
                stop(paste0("You need to manually install package '", pkg, "' from CRAN."))
            }
        }
    }
}
