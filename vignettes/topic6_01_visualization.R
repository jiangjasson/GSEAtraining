## -----------------------------------------------------------------------------
library(GSEAtraining)

lt = readRDS(system.file("extdata", "ora.rds", package = "GSEAtraining"))
diff_gene = lt$diff_gene
diff_gene = convert_to_entrez_id(diff_gene)

## -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
res = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db, pvalueCutoff = 1)
res = add_more_columns(res)

## -----------------------------------------------------------------------------
tb = res@result
head(tb)

## -----------------------------------------------------------------------------
par(mar = c(4.1, 20, 4, 1))
barplot(tb$n_hits[1:10], horiz = TRUE, names.arg = tb$Description[1:10], las = 1)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(tb[1:10, ], aes(x = n_hits, y = Description)) + 
    geom_col()

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = Description)) + 
    geom_col()

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = Description)) + 
    geom_col() +
    geom_text(aes(x = log2_fold_enrichment, label = sprintf('%1.e', p.adjust)))

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = factor(Description, levels = Description))) + 
    geom_col()

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = Description)) + 
    geom_point() +
    xlim(0, max(tb$log2_fold_enrichment[1:10]))

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = Description, size = n_hits)) + 
    geom_point() +
    xlim(0, max(tb$log2_fold_enrichment[1:10]))

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = log2_fold_enrichment, y = Description, size = n_hits, col = n_hits/gs_size)) + 
    geom_point()


## -----------------------------------------------------------------------------
ggplot(tb, aes(x = log2_fold_enrichment, y = -log10(p.adjust))) +
    geom_point(col = ifelse(tb$log2_fold_enrichment > 1 & tb$p.adjust < 0.001, "black", "grey")) +
    geom_hline(yintercept = 3, lty = 2) + geom_vline(xintercept = 1, lty = 2)

## -----------------------------------------------------------------------------
ggplot(tb, aes(x = log2_fold_enrichment, y = -log10(p.adjust), color = gs_size, size = n_hits)) +
    geom_point() + scale_colour_distiller(palette = "Spectral")


## ----fig.height = 10----------------------------------------------------------
library(enrichplot)
cnetplot(res)

## -----------------------------------------------------------------------------
heatplot(res, showCategory = 10)

## ----fig.height = 9-----------------------------------------------------------
res = pairwise_termsim(res)
treeplot(res)

## ----fig.width = 10, fig.height = 10------------------------------------------
emapplot(res)

## ----fig.height = 6, fig.width = 10-------------------------------------------
upsetplot(res)

## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "lt_enrichment_tables.rds", package = "GSEAtraining"))

set.seed(666)
terms = sample(lt[[1]]$ID, 10)

tb = NULL
for(nm in names(lt)) {
    x = lt[[nm]]
    x = x[x$ID %in% terms, colnames(x) != "geneID"]
    x$sample = nm
    tb = rbind(tb, x)
}

## -----------------------------------------------------------------------------
ggplot(tb, aes(x = sample, y = Description, size = -log10(p.adjust))) +
    geom_point(color = ifelse(tb$p.adjust < 0.05, "black", "grey"))

## -----------------------------------------------------------------------------
ggplot(tb, aes(x = sample, y = wrap_text(Description), size = -log10(p.adjust))) +
    geom_point(color = ifelse(tb$p.adjust < 0.05, "black", "grey"))

## -----------------------------------------------------------------------------
gene_diff_score = readRDS(system.file("extdata", "gene_diff_score.rds", package = "GSEAtraining"))
gene_diff_score = convert_to_entrez_id(gene_diff_score)
gene_diff_score = sort(gene_diff_score, decreasing = TRUE)
res_gsea = gseGO(geneList = gene_diff_score, ont = "BP", OrgDb = org.Hs.eg.db, pvalueCutoff = 1)
tb = res_gsea@result
head(tb)

## -----------------------------------------------------------------------------
ggplot(tb[1:10, ], aes(x = NES, y = Description)) + 
    geom_col()
ggplot(tb[1:10, ], aes(x = NES, y = Description)) + 
    geom_col(fill = ifelse(tb$NES[1:10] > 0, "red", "darkgreen"))

## -----------------------------------------------------------------------------
ggplot(tb, aes(x = NES, y = -log10(p.adjust))) +
    geom_point(col = ifelse(abs(tb$NES) > 1 & tb$p.adjust < 0.05, "black", "grey")) +
    geom_hline(yintercept = -log10(0.05), lty = 2) + geom_vline(xintercept = c(1, -1), lty = 2)

## -----------------------------------------------------------------------------
ggplot(tb, aes(x = NES, y = -log10(p.adjust), color = setSize)) +
    geom_point() + scale_colour_distiller(palette = "Spectral")

## ----fig.height = 10----------------------------------------------------------
ridgeplot(res_gsea)

## -----------------------------------------------------------------------------
gseaplot(res_gsea, geneSetID = 4)
gseaplot2(res_gsea, geneSetID = 4)

## -----------------------------------------------------------------------------
ind = c(which.max(res_gsea@result$NES), which.min(res_gsea@result$NES))
gseaplot2(res_gsea, geneSetID = ind)

