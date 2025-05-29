

.split_mat = function(mat, condition) {
    le = levels(condition)
    l_group1 = condition == le[1]
    l_group2 = !l_group1
    
    mat1 = mat[, l_group1, drop = FALSE]
    mat2 = mat[, l_group2, drop = FALSE]

    list(mat1 = mat1, mat2 = mat2)
}

#' GSEA framework
#' 
#' @param mat The input matrix.
#' @param condition The condition labels of samples. The value should be a factor where `level[1]` corresponds to the treatment group and `level[2]` corresponds to the control group.
#' 
#' @export
#' @rdname gsea_framework
#' @details
#' There are the following gene-level statistics
#' 
#' - log2fc
#' - signal to noise ratio
#' - t-value
#' - SAM t-value
#' - -log(p)
#' 
gene_level_log2fc = function(mat, condition) {
    lt = .split_mat(mat, condition)
    log2(rowMeans(lt$mat1)/rowMeans(lt$mat2))
}

#' @export
#' @rdname gsea_framework
gene_level_s2n = function(mat, condition) {
    lt = .split_mat(mat, condition)
    (rowMeans(lt$mat1) - rowMeans(lt$mat2))/(rowSds(lt$mat1) + rowSds(lt$mat2))
}

#' @export
#' @rdname gsea_framework
#' @importFrom matrixStats rowVars
gene_level_tvalue = function(mat, condition) {
    lt = .split_mat(mat, condition)
    n1 = ncol(m1)
    n2 = ncol(m2)
    v1 = rowVars(m1)
    v2 = rowVars(m2)
    sp = sqrt( ((n1-1)*v1 + (n2-1)*v2)/(n1 + n2 - 2) )*sqrt(1/n1 + 1/n2)
    (rowMeans(m1) - rowMeans(m2))/sp
}

#' @export
#' @rdname gsea_framework
#' @importFrom stats quantile
gene_level_sam = function(mat, condition) {
    lt = .split_mat(mat, condition)
    n1 = ncol(m1)
    n2 = ncol(m2)
    v1 = rowVars(m1)
    v2 = rowVars(m2)
    sp = sqrt( ((n1-1)*v1 + (n2-1)*v2)/(n1 + n2 - 2) )*sqrt(1/n1 + 1/n2)
    (rowMeans(m1) - rowMeans(m2))/(sp + quantile(sp, 0.1))
}

#' @export
#' @rdname gsea_framework
#' @importFrom genefilter rowttests
gene_level_log_p = function(mat, condition) {
    tdf = rowttests(mat, condition)  # level[1] - level[2]
    -log(tdf$p.value)*sign(tdf$statistic)
}

#' @export
#' @rdname gsea_framework
#' @param s A vector of gene scores.
#' @param l A logical vector representing whether genes in `s` are in a specific gene set.
#' 
#' @details
#' There are the following gene set level statistics:
#' 
#' - mean
#' - median
#' - maxmean
#' - wilcox ranksum 
#' - ks. Since we only support one-side test, the ES is defined as `max(abs(f2 - f1))`.
#' - chi-square. Only works when gene-level scores are binary.
set_level_mean = function(s, l) {
	mean(s[l], na.rm = TRUE)
}

#' @export
#' @rdname gsea_framework
#' @importFrom stats median
set_level_median = function(s, l) {
    median(s[l], na.rm = TRUE)
}

#' @export
#' @rdname gsea_framework
set_level_maxmean = function(s, l) {
    s = s[l]
    s1 = mean(s[s > 0], na.rm = TRUE)
    s2 = mean(s[s < 0], na.rm = TRUE)
    ifelse(s1 > abs(s2), s1, s2)
}

#' @export
#' @rdname gsea_framework
set_level_wilcox = function(s, l) {

    wilcox_stat = function(x1, x2) {
        if(length(x1) > 100) {
        x1 = sample(x1, 100)
    }
    if(length(x2) > 100) {
        x2 = sample(x2, 100)
    }
        sum(outer(x1, x2, ">"))
    }

    wilcox_stat(s[l], s[!l])
}

#' @export
#' @rdname gsea_framework
set_chi_square = function(s, l) {
    # x1: a logical vector or a binary vector
    # x2: a logical vector or a binary vector
    chisq_stat = function(x1, x2) {
        n11 = sum(x1 & x2)
        n10 = sum(x1)
        n20 = sum(!x1)
        n01 = sum(x2)
        n02 = sum(!x2)
        n = length(x1)

        n12 = n10 - n11
        n21 = n01 - n11
        n22 = n20 - n21

        p10 = n10/n
        p20 = n20/n
        p01 = n01/n
        p02 = n02/n

        e11 = n*p10*p01
        e12 = n*p10*p02
        e21 = n*p20*p01
        e22 = n*p20*p02

        stat = (n11 - e11)^2/e11 +
               (n12 - e12)^2/e12 +
               (n21 - e21)^2/e21 +
               (n22 - e22)^2/e22
        return(stat)
    }

    chisq_stat(s, l)
}

#' @export
#' @rdname gsea_framework
set_level_ks = function(s, l) {

	od = order(s, decreasing = TRUE)
	s = s[od]
	l = l[od]

	n = length(s)

    s_set = numeric(n)
    s_set[l] = abs(s[l])

    l_other = !l

    f1 = cumsum(s_set)/sum(s_set)
    f2 = cumsum(l_other)/sum(l_other)

    max(abs(f2 - f1))
}

