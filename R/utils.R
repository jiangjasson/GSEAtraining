

#' Convert gene sets between list and a data frame
#' 
#' @param lt A list of genes.
#' @rdname list_to_dataframe
#' @export
list_to_data_frame = function(lt) {
    df = data.frame(gene_set = rep(names(lt), times = sapply(lt, length)),
               gene = unlist(lt))
    rownames(df) = NULL
    df
}


#' @param df A two-column data frame.
#' 
#' @details
#' In `data_frame_to_list()`, which column contains genes and which column contains gene sets
#' are automatically checked by the number of genes and gene sets. Basically number of genes
#' should be larger than the number of gene sets.
#' 
#' @rdname list_to_dataframe
#' @export
data_frame_to_list = function(df) {
    n1 = length(unique(df[, 1]))
    n2 = length(unique(df[, 2]))
    if(n1 < n2) {
        split(df[, 2], df[, 1])
    } else {
        split(df[, 1], df[, 2])
    }
}



#' Wrap long text into multiple lines
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



# See \url{https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes,_unequal_variances_(sX1_%3E_2sX2_or_sX2_%3E_2sX1)}.
effective_degrees_of_freedom = function(m1, m2) {
    v1 = rowVars(m1)
    v2 = rowVars(m2)

    n1 = ncol(m1)
    n2 = ncol(m2)

    ( v1/n1 + v2/n2 )^2/( (v1/n1)^2/(n1 - 1) + (v2/n2)^2/(n2 - 1) )
}


#' Create a GRanges object from seqlengths
#' 
#' @param x Any object that supports the `seqlengths()` method.
#' @param seqnames A vector seqnames.
#' 
#' @details
#' It basically generates a `GRanges` treating whole chromosomes as genomic intervals.
#' 
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' seqlengths_as_GRanges(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' seqlengths_as_GRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'     paste0("chr", c(1:22, "X", "Y")))
seqlengths_as_GRanges = function(x, seqnames = NULL) {
    len = seqlengths(x)
    if(!is.null(seqnames)) {
        len = len[names(len) %in% as.vector(seqnames)]
    }
    GRanges(seqnames = names(len), ranges = IRanges(1, len))
}
