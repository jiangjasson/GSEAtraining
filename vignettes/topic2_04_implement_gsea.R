## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "p53_expr.rds", package = "GSEAtraining"))
expr = lt$expr
condition = lt$condition

ln = strsplit(readLines(system.file("extdata", "c2.symbols.gmt", package = "GSEAtraining")), "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])

geneset = gs[["p53hypoxiaPathway"]]

## -----------------------------------------------------------------------------
length(geneset)

## -----------------------------------------------------------------------------
s = apply(expr, 1, function(x) {
	x1 = x[condition == "WT"]
	x2 = x[condition == "MUT"]
	(mean(x1) - mean(x2))/(sd(x1) + sd(x2))
})

## -----------------------------------------------------------------------------
s = sort(s, decreasing = TRUE)

## -----------------------------------------------------------------------------
## original GSEA
l_set = names(s) %in% geneset
f1 = cumsum(l_set)/sum(l_set)

## or
binary_set = l_set + 0
f1 = cumsum(binary_set)/sum(binary_set)

## -----------------------------------------------------------------------------
l_other = ! names(s) %in% geneset
f2 = cumsum(l_other)/sum(l_other)

## -----------------------------------------------------------------------------
n = length(s)
plot(1:n, f1, type = "l", col = "red", xlab = "sorted genes")
lines(1:n, f2, col = "blue")

## -----------------------------------------------------------------------------
plot(f1 - f2, type = "l", xlab = "sorted genes")
abline(h = 0, lty = 2, col = "grey")
points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")

## -----------------------------------------------------------------------------
es = max(f1 - f2)
es

## -----------------------------------------------------------------------------
plot(f1 - f2, type = "l", xlab = "sorted genes")
abline(h = 0, lty = 2, col = "grey")
points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")
abline(v = which.max(f1 - f2), lty = 3, col = "blue")

## -----------------------------------------------------------------------------
ks.test(which(l_set), which(l_other))

## -----------------------------------------------------------------------------
library(matrixStats)
# expr: the complete expression matrix
# condition: the condition labels of samples
# cmp: a vector of two, cmp[1] - cmp[2] > 0 means up-regulation
# geneset: A vector of genes
calculate_es = function(expr, condition, cmp, geneset) {

	m1 = expr[, condition == cmp[1]]  # only samples in group 1
	m2 = expr[, condition == cmp[2]]  # only samples in group 2

	s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2)) # a gene-level difference socre (S2N ratio) 

	s = sort(s, decreasing = TRUE)  # ranked gene list

	l_set = names(s) %in% geneset
	f1 = cumsum(l_set)/sum(l_set)   # CDF for genes in the set

	l_other = !l_set
	f2 = cumsum(l_other)/sum(l_other)  # CDF for genes not in the set

	max(f1 - f2)
}

## -----------------------------------------------------------------------------
es = calculate_es(expr, condition, cmp = c("WT", "MUT"), geneset = geneset)
es

## -----------------------------------------------------------------------------
set.seed(123)
es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es(expr, sample(condition), 
	    cmp = c("WT", "MUT"), geneset = geneset)
}

## -----------------------------------------------------------------------------
sum(es_rand >= es)/1000

## -----------------------------------------------------------------------------
hist(es_rand)
abline(v = es, col = "red")

## -----------------------------------------------------------------------------
calculate_es_v2 = function(expr, condition, cmp, geneset, plot = FALSE, power = 1) {

    m1 = expr[, condition == cmp[1]]
    m2 = expr[, condition == cmp[2]]

    s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

    s = sort(s, decreasing = TRUE)

    l_set = names(s) %in% geneset
    
    # f1 = cumsum(l_set)/sum(l_set)  # <<-- the original line
    
    s_set = abs(s)^power   # <<-- here
    s_set[!l_set] = 0
    f1 = cumsum(s_set)/sum(s_set)  ## <<- here

    l_other = !l_set
    f2 = cumsum(l_other)/sum(l_other)

    if(plot) {
        plot(f1 - f2, type = "l", xlab = "sorted genes")
        abline(h = 0, lty = 2, col = "grey")
        points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")
        abline(v = which.max(f1 - f2), lty = 3, col = "blue")
    }

    max(f1 - f2)
}

## -----------------------------------------------------------------------------
es = calculate_es_v2(expr, condition, cmp = c("WT", "MUT"), plot = TRUE, 
    geneset = geneset)

## ----results = "none", fig.width = 10-----------------------------------------
par(mfrow = c(1, 2))
calculate_es_v2(expr, condition, cmp = c("WT", "MUT"), plot = TRUE, power = 0, 
    geneset = geneset)  # same as the original GSEA
title("power = 0")
calculate_es_v2(expr, condition, cmp = c("WT", "MUT"), plot = TRUE, power = 2, 
    geneset = geneset)
title("power = 2")

## -----------------------------------------------------------------------------
es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es_v2(expr, sample(condition), 
	    cmp = c("WT", "MUT"), geneset = geneset)
}

## -----------------------------------------------------------------------------
sum(es_rand >= es)/1000

## -----------------------------------------------------------------------------
hist(es_rand, xlim = c(0, 1))
abline(v = es, col = "red")

## -----------------------------------------------------------------------------
plot(abs(s))

## -----------------------------------------------------------------------------
# s: a vector of pre-calcualted gene-level scores
# s should be sorted
calculate_es_v2_gene_perm = function(s, geneset, perm = FALSE, plot = FALSE, power = 1) {
	
	if(perm) {
	    # s is still sorted, but the gene labels are randomly shuffled
	    # to break the associations between gene scores and gene labels.
		names(s) = sample(names(s))  ## <<- here
	}

	l_set = names(s) %in% geneset
	s_set = abs(s)^power
	s_set[!l_set] = 0
	f1 = cumsum(s_set)/sum(s_set)

	l_other = !l_set
	f2 = cumsum(l_other)/sum(l_other)

	if(plot) {
        plot(f1 - f2, type = "l", xlab = "sorted genes")
        abline(h = 0, lty = 2, col = "grey")
        points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")
        abline(v = which.max(f1 - f2), lty = 3, col = "blue")
    }

	max(f1 - f2)
}

## -----------------------------------------------------------------------------
# pre-calculate gene-level scores
m1 = expr[, condition == "WT"]
m2 = expr[, condition == "MUT"]

s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))
s = sort(s, decreasing = TRUE)  # must be pre-sorted

## -----------------------------------------------------------------------------
es = calculate_es_v2_gene_perm(s, geneset, plot = TRUE)

## -----------------------------------------------------------------------------
es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es_v2_gene_perm(s, geneset, perm = TRUE)
}

sum(es_rand >= es)/1000

## -----------------------------------------------------------------------------
hist(es_rand, xlim = c(0, 1))
abline(v = es, col = "red")

## -----------------------------------------------------------------------------
es/mean(es_rand)

## -----------------------------------------------------------------------------
(es - mean(es_rand))/sd(es_rand)

## -----------------------------------------------------------------------------
library(GSEAtraining)
gs_hallmark = get_msigdb(version = "2023.2.Hs", collection = "h.all", gene_id_type = "symbol")
n_gs = length(gs_hallmark)

## -----------------------------------------------------------------------------
p1 = numeric(n_gs)
es1 = numeric(n_gs)
nes1 = numeric(n_gs)
z1 = numeric(n_gs)

for(i in 1:n_gs) {
	print(i)
	geneset = gs_hallmark[[i]]

	es1[i] = calculate_es_v2(expr, condition, 
		    cmp = c("WT", "MUT"), geneset = geneset)

	es_rand = numeric(1000)
	for(k in 1:1000) {
		es_rand[k] = calculate_es_v2(expr, sample(condition), 
		    cmp = c("WT", "MUT"), geneset = geneset)
	}
	p1[i] = sum(es_rand >= es1[i])/1000
	nes1[i] = es1[i]/mean(es_rand)
	z1[i] = (es1[i] - mean(es_rand))/sd(es_rand)
}

df1 = data.frame(
	geneset = names(gs_hallmark),
	es = es1,
	nes = nes1,
	z = z1,
	p_value = p1
)

## -----------------------------------------------------------------------------
m1 = expr[, condition == "WT"]
m2 = expr[, condition == "MUT"]

s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))
s = sort(s, decreasing = TRUE)  # must be pre-sorted

p2 = numeric(n_gs)
es2 = numeric(n_gs)
nes2 = numeric(n_gs)
z2 = numeric(n_gs)

for(i in 1:n_gs) {
	print(i)
	geneset = gs_hallmark[[i]]

	es2[i] = calculate_es_v2_gene_perm(s, geneset = geneset)
	
	es_rand = numeric(1000)
	for(k in 1:1000) {
		es_rand[k] = calculate_es_v2_gene_perm(s, geneset = geneset, perm = TRUE)
	}
	p2[i] = sum(es_rand >= es2[i])/1000
	nes2[i] = es2[i]/mean(es_rand)
	z2[i] = (es2[i] - mean(es_rand))/sd(es_rand)
}

df2 = data.frame(
	geneset = names(gs_hallmark),
	es = es2,
	nes = nes2,
	z = z2,
	p_value = p2
)

## ----fig.width = 10, fig.height = 10------------------------------------------
par(mfrow = c(2, 2))
plot(df1$es, df2$es, main = "ES",
	xlab = "sample permutation", ylab = "gene permutation")
plot(df1$nes, df2$nes, main = "NES", xlim = c(0, 2.5), ylim = c(0, 2.5),
	xlab = "sample permutation", ylab = "gene permutation")
plot(df1$z, df2$z, main = "Z-score",
	xlab = "sample permutation", ylab = "gene permutation")
plot(df1$p_value, df2$p_value, main = "P-value",
	xlab = "sample permutation", ylab = "gene permutation")

