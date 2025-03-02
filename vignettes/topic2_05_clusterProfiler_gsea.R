## -----------------------------------------------------------------------------
gene_diff_score = readRDS(system.file("extdata", "gene_diff_score.rds", package = "GSEAtraining"))
head(gene_diff_score)

## -----------------------------------------------------------------------------
library(org.Hs.eg.db)
map = mapIds(org.Hs.eg.db, keys = names(gene_diff_score), keytype = "SYMBOL", column = "ENTREZID")
head(map)

# names(gene_diff_score) = map[ names(gene_diff_score) ]

## -----------------------------------------------------------------------------
library(GSEAtraining)
gene_diff_score = convert_to_entrez_id(gene_diff_score)
head(gene_diff_score)

## ----message=FALSE------------------------------------------------------------
library(clusterProfiler)

## -----------------------------------------------------------------------------
gene_diff_score = sort(gene_diff_score, decreasing = TRUE)

## ----message = FALSE----------------------------------------------------------
library(org.Hs.eg.db)
tb = gseGO(geneList = gene_diff_score, ont = "BP", OrgDb = org.Hs.eg.db)
head(tb)

## ----message = FALSE----------------------------------------------------------
tb = gseKEGG(geneList = gene_diff_score, organism = "hsa")
head(tb)

## -----------------------------------------------------------------------------
library(ReactomePA)
tb = gsePathway(geneList = gene_diff_score, organism = "human")
head(tb)

## -----------------------------------------------------------------------------
library(DOSE)
tb = gseDO(geneList = gene_diff_score)
head(tb)

## -----------------------------------------------------------------------------
library(msigdbr)
gene_sets = msigdbr(category = "H")
map = gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene = as.character(map$entrez_gene)

tb = GSEA(geneList = gene_diff_score, TERM2GENE = map)
head(tb)

## -----------------------------------------------------------------------------
library(org.Ss.eg.db)
all_genes = keys(org.Ss.eg.db, keytype = "ENTREZID")
scores = sort(rnorm(length(all_genes)), decreasing = TRUE)
names(scores) = sample(all_genes)

## ----message = FALSE----------------------------------------------------------
tb = gseGO(geneList = scores, ont = "BP", OrgDb = org.Ss.eg.db, pvalueCutoff = 1)
head(tb)

## ----message = FALSE----------------------------------------------------------
tb = gseKEGG(geneList = scores, organism = "ssc", pvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
gene_sets = msigdbr(species = "pig", category = "H")
map = gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene = as.character(map$entrez_gene)

tb = GSEA(geneList = scores, TERM2GENE = map, pvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()
org_db = ah[["AH118180"]]

all_genes = keys(org_db, keytype = "ENTREZID")
scores = sort(rnorm(length(all_genes)), decreasing = TRUE)
names(scores) = sample(all_genes)

tb = gseGO(geneList = scores, ont = "BP", OrgDb = org_db, pvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
class(tb)
str(tb)

## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "p53_expr.rds", package = "GSEAtraining"))
expr = lt$expr
condition = lt$condition

## -----------------------------------------------------------------------------
mat = convert_to_entrez_id(expr)

## -----------------------------------------------------------------------------
ind_MUT = which(condition == "MUT")
ind_WT = which(condition == "WT")

## -----------------------------------------------------------------------------
s2n = apply(mat, 1, function(x) {
    x1 = x[ind_MUT]
    x2 = x[ind_WT]

    (mean(x1) - mean(x2))/(sd(x1) + sd(x2))
})

tvalue = apply(mat, 1, function(x) {
    x1 = x[ind_MUT]
    x2 = x[ind_WT]

    t.test(x1, x2)$statistic
})

log2fc = apply(mat, 1, function(x) {
    x1 = x[ind_MUT]
    x2 = x[ind_WT]

    log2(mean(x1)/mean(x2))
})

logp = apply(mat, 1, function(x) {
    x1 = x[ind_MUT]
    x2 = x[ind_WT]

    t = t.test(x1, x2)
    sign(t$statistic)*(-log10(t$p.value))
})

## -----------------------------------------------------------------------------
s2n = sort(s2n, decreasing = TRUE)
tvalue = sort(tvalue, decreasing = TRUE)
log2fc = sort(log2fc, decreasing = TRUE)
logp = sort(logp, decreasing = TRUE)

## -----------------------------------------------------------------------------
tb1 = gseGO(geneList = s2n, ont = "BP", OrgDb = org.Hs.eg.db)
tb2 = gseGO(geneList = tvalue, ont = "BP", OrgDb = org.Hs.eg.db)
tb3 = gseGO(geneList = log2fc, ont = "BP", OrgDb = org.Hs.eg.db)
tb4 = gseGO(geneList = logp, ont = "BP", OrgDb = org.Hs.eg.db)

## -----------------------------------------------------------------------------
tb1 = tb1@result
tb2 = tb2@result
tb3 = tb3@result
tb4 = tb4@result

## -----------------------------------------------------------------------------
library(eulerr)
plot(euler(list("s2n" = tb1$ID[tb1$p.adjust < 0.05],
                "tvalue" = tb2$ID[tb2$p.adjust < 0.05],
                "log2fc" = tb3$ID[tb3$p.adjust < 0.05],
                "logp" = tb4$ID[tb4$p.adjust < 0.05])),
    quantities = TRUE)

## ----fig.width = 8, fig.height = 8--------------------------------------------
par(mfrow = c(2, 2))
plot(s2n)
plot(tvalue)
plot(log2fc)
plot(logp)

