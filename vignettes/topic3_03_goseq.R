## -----------------------------------------------------------------------------
tb = readRDS(system.file("extdata", "de.rds", package = "GSEAtraining"))
tb = tb[tb$symbol != "", ]  # keep pc genes
tb = tb[!is.na(tb$p_value), ]
tb$fdr = p.adjust(tb$p_value, "BH")

genes = ifelse(tb$fdr < 0.05, 1, 0)
names(genes) = tb$ensembl
head(genes)
table(genes)

## ----message = FALSE----------------------------------------------------------
library(goseq)
pwf = nullp(genes, "hg19", "ensGene")

## ----message = FALSE----------------------------------------------------------
tb1 = goseq(pwf, "hg19", "ensGene")
head(tb1)

## -----------------------------------------------------------------------------
tb = readRDS(system.file("extdata", "de.rds", package = "GSEAtraining"))
tb = tb[tb$symbol != "", ]  # keep pc genes

tb$fdr = p.adjust(tb$p_value, "BH")
genes = ifelse(tb$fdr < 0.05, 1, 0)
genes[is.na(genes)] = 0
names(genes) = tb$ensembl
table(genes)

## -----------------------------------------------------------------------------
pwf = nullp(genes, "hg19", "ensGene")

## ----message = FALSE----------------------------------------------------------
tb2 = goseq(pwf, "hg19", "ensGene")

## -----------------------------------------------------------------------------
library(eulerr)
lt = list(
    test1 = tb1$category[p.adjust(tb1$over_represented_pvalue, "BH") < 0.05], 
    test2 = tb2$category[p.adjust(tb2$over_represented_pvalue, "BH") < 0.05]
)
plot(euler(lt), quantities = TRUE)

## -----------------------------------------------------------------------------
l = tb$fdr < 0.05; l[is.na(l)] = FALSE
diff_gene = tb$ensembl[l]

## -----------------------------------------------------------------------------
library(GSEAtraining)
diff_gene = convert_to_entrez_id(diff_gene)

## ----message = FALSE----------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
tb3 = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db)

## -----------------------------------------------------------------------------
tb1 = tb1[tb1$ontology == "BP", ]
tb2 = tb2[tb2$ontology == "BP", ]

library(eulerr)
lt = list(
    test1 = tb1$category[p.adjust(tb1$over_represented_pvalue, "BH") < 0.05], 
    test2 = tb2$category[p.adjust(tb2$over_represented_pvalue, "BH") < 0.05],
    ora = tb3$ID[p.adjust(tb3$pvalue, "BH") < 0.05]
)
plot(euler(lt), quantities = TRUE)

