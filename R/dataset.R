
#' P53 dataset
#' 
#' @usage data(p53_dataset)
#' 
#' @details
#' P53 dataset as well as the C2 gene sets are from \url{https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/}.
#' 
#' The following code is used to generate the dataset:
#' 
#' ```
#' expr = read.table(url("https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/P53_collapsed_symbols.gct"), 
#'     skip = 2, header = TRUE, row.names = 1, sep = "\t", quote = "")
#' expr = as.matrix(expr[, -1])
#' 
#' condition = readLines(url("https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/P53.cls"))[3]
#' condition = strsplit(condition, " ")[[1]]
#' condition = factor(condition, levels = c("MUT", "WT"))
#' 
#' ln = readLines(url("https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/c2.symbols.gmt"))
#' lt = strsplit(ln, "\t")
#' 
#' gs = lapply(lt, function(x) x[-(1:2)])
#' names(gs) = sapply(lt, function(x) x[1])
#' 
#' library(matrixStats)
#' m1 = expr[, condition == "MUT"] 
#' m2 = expr[, condition == "WT"]  
#' s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2)) 
#' s = sort(s, decreasing = TRUE)
#' 
#' p53_dataset = list(expr = expr, condition = condition, s2n = s, gs = gs)
#' ```
#' 
#' A positive signal-to-noise ratio (`s2n`) means up-regulation in `MUT`.
#' 
"p53_dataset"
