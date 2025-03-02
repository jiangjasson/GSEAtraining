## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "p53_expr.rds", package = "GSEAtraining"))
expr = lt$expr
condition = lt$condition

ln = strsplit(readLines(system.file("extdata", "c2.symbols.gmt", package = "GSEAtraining")), "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])

## -----------------------------------------------------------------------------
library(GSVA)
gs_mat = gsva(gsvaParam(expr, gs), verbose = FALSE)
dim(gs_mat)

## ----message = FALSE----------------------------------------------------------
library(ComplexHeatmap)
Heatmap(gs_mat, top_annotation = HeatmapAnnotation(cond = condition),
    column_split = condition)

## ----message=FALSE------------------------------------------------------------
library(genefilter)
tdf = rowttests(gs_mat, factor(condition))
tdf$fdr = p.adjust(tdf$p.value, "BH")

## -----------------------------------------------------------------------------
tdf[tdf$fdr < 0.05, ]

## -----------------------------------------------------------------------------
s = rowttests(expr, factor(condition))[, "statistic"]
names(s) = rownames(expr)

## ----message = FALSE----------------------------------------------------------
map = data.frame(
    gene_set = rep(names(gs), times = sapply(gs, length)),
    gene = unlist(gs)
)
    
library(clusterProfiler)
tb = GSEA(geneList = sort(s, decreasing = TRUE), TERM2GENE = map, pvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
g = intersect(gs[["p53Pathway"]], rownames(expr))
mm = expr[g, ]
mm = t(scale(t(mm)))  # z-score transformation

Heatmap(mm, name = "z-score", top_annotation = HeatmapAnnotation(cond = condition),
    column_title = "p53Pathway",
    right_annotation = rowAnnotation(bar = anno_barplot(s[g], axis_param = list(direction = "reverse"), width = unit(2, "cm"))),
    column_split = condition) %v%
Heatmap(gs_mat["hivnefPathway", , drop = FALSE], name = "set variation")

## ----fig.width = 10, fig.height = 5-------------------------------------------
cn = intersect(rownames(tdf), tb$ID)
par(mfrow = c(1, 2))
plot(tdf[cn, "p.value"], tb[cn, "pvalue"],
    xlab = "GSVA", ylab = "normal GSEA", main = "p-values")
plot(tdf[cn, "statistic"], tb[cn, "NES"],
    xlab = "GSVA / t-values", ylab = "normal GSEA / NES")
abline(h = 0, lty = 2, col = "grey")
abline(v = 0, lty = 2, col = "grey")

## -----------------------------------------------------------------------------
g = intersect(gs[["hivnefPathway"]], rownames(expr))
mm = expr[g, ]
mm = t(scale(t(mm)))  # z-score transformation

Heatmap(mm, name = "z-score", top_annotation = HeatmapAnnotation(cond = condition),
    right_annotation = rowAnnotation(bar = anno_barplot(s[g], axis_param = list(direction = "reverse"), width = unit(2, "cm"))),
    column_split = condition) %v%
Heatmap(gs_mat["hivnefPathway", , drop = FALSE], name = "set variation")

