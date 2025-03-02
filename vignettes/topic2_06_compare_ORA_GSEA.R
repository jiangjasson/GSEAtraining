## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "p53_expr.rds", package = "GSEAtraining"))
expr = lt$expr
condition = lt$condition

## -----------------------------------------------------------------------------
library(GSEAtraining)
expr = convert_to_entrez_id(expr)

## -----------------------------------------------------------------------------
p = apply(expr, 1, function(x) {
    x1 = x[condition == "WT"]
    x2 = x[condition == "MUT"]
    
    t.test(x1, x2)$p.value
})

## -----------------------------------------------------------------------------
library(genefilter)
tdf = rowttests(expr, factor(condition))  # the second must be a "factor"
tdf$fdr = p.adjust(tdf$p.value, "BH")
sum(tdf$fdr < 0.05)  # number of diff genes

## -----------------------------------------------------------------------------
plot(sort(tdf$statistic))

## -----------------------------------------------------------------------------
sum(abs(tdf$statistic) > 2)

## ----message = FALSE----------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
l_sig = abs(tdf$statistic) > 2
tb_ora = enrichGO(gene  = rownames(expr)[l_sig],
                  OrgDb = org.Hs.eg.db,
                  ont   = "BP",
                  pAdjustMethod = "BH")
tb_ora = as.data.frame(tb_ora)  # only significant ones

## -----------------------------------------------------------------------------
s = tdf$statistic
names(s) = rownames(tdf)   # s must have names (gene IDs)
s = sort(s, decreasing = TRUE)  # s must be pre-sorted
tb_gsea = gseGO(geneList = s, 
                OrgDb = org.Hs.eg.db,
                ont   = "BP",
                pAdjustMethod = "BH")
tb_gsea = as.data.frame(tb_gsea)  # only significant ones

## -----------------------------------------------------------------------------
library(eulerr)
plot(euler(list(ORA = tb_ora$ID, GSEA = tb_gsea$ID)), quantities = TRUE)

