## -----------------------------------------------------------------------------
library(org.Hs.eg.db)

## -----------------------------------------------------------------------------
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

## -----------------------------------------------------------------------------
genes = c("TP53", "MDM2")

## -----------------------------------------------------------------------------
map = select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = "ENTREZID")
map

## -----------------------------------------------------------------------------
map = select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = "ENSEMBL")
map

## -----------------------------------------------------------------------------
select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))

## -----------------------------------------------------------------------------
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID")
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENSEMBL")

## -----------------------------------------------------------------------------
genes = c("TP53", "MMD2")
map = select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = "ENTREZID")
map

## -----------------------------------------------------------------------------
map = select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("ENTREZID", "GENETYPE"))
map

## -----------------------------------------------------------------------------
map[map$GENETYPE == "protein-coding", ]

## -----------------------------------------------------------------------------
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "asNA")
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "list")
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "filter")

gene_type = as.list(org.Hs.egGENETYPE)
# this is very slow
mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", 
    multiVals = function(x) {
        x2 = x[sapply(as.list(gene_type[x]), function(y) y == "protein-coding")]
        if(length(x2)) {
            x2[1]
        } else {
            NA
        }
})

## -----------------------------------------------------------------------------
ls(envir = asNamespace("org.Hs.eg.db"))

## -----------------------------------------------------------------------------
org.Hs.egSYMBOL

## -----------------------------------------------------------------------------
org.Hs.egSYMBOL2EG[["TP53"]]

## -----------------------------------------------------------------------------
lt = as.list(org.Hs.egSYMBOL2EG)
lt[genes]

# or
as.list(org.Hs.egSYMBOL2EG[genes])

## -----------------------------------------------------------------------------
tb = toTable(org.Hs.egSYMBOL2EG)

## ----eval = FALSE-------------------------------------------------------------
# library(GSEAtraining)
# convert_to_entrez_id(genes)

## -----------------------------------------------------------------------------
lt = readRDS(system.file("extdata", "ora.rds", package = "GSEAtraining"))
diff_gene = lt$diff_gene
head(diff_gene)

## ----eval = FALSE-------------------------------------------------------------
# map = select(org.Hs.eg.db, keys = diff_gene, keytype = "SYMBOL", columns = "ENTREZID")
# unique(map$ENTREZID)
# 
# # or
# map = mapIds(org.Hs.eg.db, keys = diff_gene, keytype = "SYMBOL", column = "ENTREZID")

## -----------------------------------------------------------------------------
all_symbols = keys(org.Hs.eg.db, keytype = "SYMBOL")

# use select
map = select(org.Hs.eg.db, keys = all_symbols, keytype = "SYMBOL", columns = "ENTREZID")

tb = table(map$SYMBOL)
sum(tb > 1)/length(tb)

# use mapIds
map = mapIds(org.Hs.eg.db, keys = all_symbols, keytype = "SYMBOL", column = "ENTREZID", 
    multiVals = "list")
n = sapply(map, length)
sum(n > 1)/length(n)

# use org.Hs.egSYMBOL2EG
lt = as.list(org.Hs.egSYMBOL2EG)
n = sapply(lt, length)
sum(n > 1)/length(n)

## -----------------------------------------------------------------------------
# use select
map = select(org.Hs.eg.db, keys = all_symbols, keytype = "SYMBOL", 
    columns = c("ENTREZID", "GENETYPE"))
map = map[map$GENETYPE == "protein-coding", ]

tb = table(map$SYMBOL)
sum(tb > 1)/length(tb)

## -----------------------------------------------------------------------------
all_ids = keys(org.Hs.eg.db, keytype = "ENSEMBL")

# use select
map = select(org.Hs.eg.db, keys = all_ids, keytype = "ENSEMBL", columns = "ENTREZID")

tb = table(map$ENSEMBL)
sum(tb > 1)/length(tb)

tail(sort(tb))

## -----------------------------------------------------------------------------
# use select
map = select(org.Hs.eg.db, keys = all_ids, keytype = "ENSEMBL", columns = c("ENTREZID", "GENETYPE"))
map = map[map$GENETYPE == "protein-coding", ]

tb = table(map$ENSEMBL)
sum(tb > 1)/length(tb)

## -----------------------------------------------------------------------------
map = select(org.Hs.eg.db, keys = all_ids, keytype = "ENSEMBL", 
    columns = c("ENTREZID", "SYMBOL", "GENETYPE"))
map[map$ENSEMBL %in% names(tail(sort(tb))), ]

