## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "ora.rds", package = "GSEAtraining"))
diff_gene = lt$diff_gene
bg_gene = lt$bg_gene
head(diff_gene)

## ----message=FALSE------------------------------------------------------------
library(clusterProfiler)

## -----------------------------------------------------------------------------
library(GSEAtraining)
diff_gene = convert_to_entrez_id(diff_gene)
head(diff_gene)
length(diff_gene)

## ----message = FALSE----------------------------------------------------------
library(org.Hs.eg.db)
tb = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db)
head(tb)

## ----message = FALSE----------------------------------------------------------
tb = enrichKEGG(gene = diff_gene, organism = "hsa")
head(tb)

## -----------------------------------------------------------------------------
library(ReactomePA)
tb = enrichPathway(gene = diff_gene, organism = "human")
head(tb)

## -----------------------------------------------------------------------------
library(DOSE)
tb = enrichDO(gene = diff_gene)
head(tb)

## -----------------------------------------------------------------------------
library(msigdbr)
gene_sets = msigdbr(category = "H")
map = gene_sets[, c("gs_name", "entrez_gene")]
head(map)

tb = enricher(gene = diff_gene, TERM2GENE = map)
head(tb)

## -----------------------------------------------------------------------------
bg_gene = convert_to_entrez_id(bg_gene)
bg_gene = sample(bg_gene, 10000)
tb = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db, universe = bg_gene)
head(tb)

## -----------------------------------------------------------------------------
tb1 = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db)
tb2 = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db, universe = bg_gene)

tb1 = tb1@result
tb2 = tb2@result
cn = intersect(tb1$ID, tb2$ID)

plot(tb1[cn, "pvalue"], tb2[cn, "pvalue"], log = "xy",
	xlim = c(1e-25, 1), ylim = c(1e-25, 1),
	xlab = "default background", ylab = "self-defined background")

## -----------------------------------------------------------------------------
library(eulerr)
plot(euler(list("default_background" = tb1$ID[tb1$p.adjust < 0.05],
	            "bg_gene" = tb2$ID[tb2$p.adjust < 0.05])),
    quantities = TRUE)

## -----------------------------------------------------------------------------
library(org.Ss.eg.db)
diff_gene = random_genes(org.Ss.eg.db, 1000, "ENTREZID")

## ----message = FALSE----------------------------------------------------------
tb = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Ss.eg.db, 
	pvalueCutoff = 1, qvalueCutoff = 1)
head(tb)

## ----message = FALSE----------------------------------------------------------
tb = enrichKEGG(gene = diff_gene, organism = "ssc", pvalueCutoff = 1, qvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
gene_sets = msigdbr(species = "pig", category = "H")
map = gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene = as.character(map$entrez_gene)

tb = enricher(gene = diff_gene, TERM2GENE = map, pvalueCutoff = 1, qvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()
org_db = ah[["AH118180"]]

diff_gene = random_genes(org_db, 1000, "ENTREZID")
tb = enrichGO(diff_gene, ont = "BP", OrgDb = org_db, 
	pvalueCutoff = 1, qvalueCutoff = 1)
head(tb)

## -----------------------------------------------------------------------------
class(tb)
str(tb)

## -----------------------------------------------------------------------------
tb2 = add_more_columns(tb)
head(tb2@result)

## -----------------------------------------------------------------------------
genes = readRDS(system.file("extdata", "webgestalt_example_gene_list.rds", 
	package = "GSEAtraining"))

