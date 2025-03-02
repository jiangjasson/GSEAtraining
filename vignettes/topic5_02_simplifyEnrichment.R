## ----echo = FALSE-------------------------------------------------------------
library(simplifyEnrichment)
se_opt$verbose = FALSE

## ----fig.width = 8, fig.height = 4.5------------------------------------------
library(simplifyEnrichment)
set.seed(888)
go_id = random_GO(500)
df = simplifyGO(go_id)
head(df)

## ----echo = FALSE, message = FALSE, fig.width = 8, fig.height = 7-------------
library(cola)
data(golub_cola) 
res = golub_cola["ATC:skmeans"]
get_signatures(res, k = 3)

## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "lt_enrichment_tables.rds", package = "GSEAtraining"))
length(lt)
head(lt[[1]])

## ----fig.width = 10, fig.height = 6, eval = FALSE-----------------------------
# simplifyGOFromMultipleLists(lt, padj_cutoff = 0.001)

## ----fig.width = 8, fig.height = 4.5, eval = FALSE----------------------------
# df = lt[[1]]
# go_id = df$ID[df$p.adjust < 0.01]
# simplifyGO(go_id)

## ----fig.width = 8, fig.height = 5, eval = FALSE------------------------------
# l = df$p.adjust < 0.01
# summarizeGO(df$ID[l], -log10(df$p.adjust)[l], axis_label = "average -log10(p.adjust)")

## ----fig.width = 8, fig.height = 5, eval = FALSE------------------------------
# l = df$p.adjust < 0.01
# summarizeGO(df$ID[l], df$log2_fold_enrichment[l], axis_label = "average log2(fold enrichment)")

## ----fig.width = 8, fig.height = 7, eval = FALSE------------------------------
# value = lapply(lt, function(df) {
#     v = -log10(df$p.adjust)
#     names(v) = df$ID
#     v[df$p.adjust < 0.001]
# })
# summarizeGO(value = value, axis_label = "average -log10(p.adjust)",
#     legend_title = "-log10(p.adjust)")

## ----fig.width = 8, fig.height = 7, eval = FALSE------------------------------
# value = lapply(lt, function(df) {
#     v = df$log2_fold_enrichment
#     names(v) = df$ID
#     v[df$p.adjust < 0.001]
# })
# summarizeGO(value = value, axis_label = "average log2_fold_enrichment",
#     legend_title = "log2(fold_enrichment)")

## -----------------------------------------------------------------------------
de = readRDS(system.file("extdata", "de.rds", package = "GSEAtraining"))
de = de[, c("symbol", "p_value")]
de = de[!is.na(de$p_value), ]

de_genes_1 = de$symbol[de$p_value < 0.05]
de_genes_2 = de$symbol[de$p_value < 0.01]
de_genes_3 = de$symbol[de$p_value < 0.001]

library(GSEAtraining)
library(org.Hs.eg.db)
gs = get_GO_gene_sets_from_orgdb(org.Hs.eg.db, "BP", gene_id_type = "SYMBOL")

tb1 = ora(de_genes_1, gs)
tb2 = ora(de_genes_2, gs)
tb3 = ora(de_genes_3, gs)

## ----eval = FALSE-------------------------------------------------------------
# lt = list(
#     "DE_0.05" = tb1,
#     "DE_0.01" = tb2,
#     "DE_0.001" = tb3
# )
# 
# simplifyGOFromMultipleLists(lt, padj_cutoff = 0.05)

