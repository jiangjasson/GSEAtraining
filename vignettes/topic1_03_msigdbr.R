## -----------------------------------------------------------------------------
download.file("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
  destfile = "h.all.v2023.2.Hs.symbols.gmt")

## -----------------------------------------------------------------------------
ln = readLines("h.all.v2023.2.Hs.symbols.gmt")
ln = strsplit(ln, "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])
gs[1:2]

## -----------------------------------------------------------------------------
library(rGREAT)
gs = read_gmt("h.all.v2023.2.Hs.symbols.gmt")

gs2 = read_gmt("h.all.v2023.2.Hs.symbols.gmt", 
    from = "SYMBOL", to = "ENSEMBL", orgdb = "org.Hs.eg.db")
gs2[1:2]

## ----echo = FALSE-------------------------------------------------------------
invisible(file.remove("h.all.v2023.2.Hs.symbols.gmt"))

## ----message = FALSE----------------------------------------------------------
library(msigdbr)

## -----------------------------------------------------------------------------
msigdbr_species()

## -----------------------------------------------------------------------------
as.data.frame(msigdbr_collections())

## -----------------------------------------------------------------------------
gene_sets = msigdbr(category = "C2", subcategory = "CP:KEGG")
gene_sets

## -----------------------------------------------------------------------------
gene_sets = as.data.frame(gene_sets)
head(gene_sets)

## -----------------------------------------------------------------------------
gene_sets[, c("gs_name", "ensembl_gene")] |> head()

## -----------------------------------------------------------------------------
gene_sets = msigdbr(species = "dog", category = "C2", subcategory = "CP:BIOCARTA")
gene_sets = as.data.frame(gene_sets)
head(gene_sets)

## -----------------------------------------------------------------------------
gene_sets = msigdbr(category = "C5", subcategory = "GO:BP")

tb = gene_sets[, c("gs_name", "gene_symbol")]
dim(tb)
dim(unique(tb))

## -----------------------------------------------------------------------------
library(GSEAtraining)
list_msigdb_versions()
list_msigdb_collections("2023.2.Hs")
lt = get_msigdb(version = "2023.2.Hs", collection = "h.all")
lt[1:2]

## -----------------------------------------------------------------------------
download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=COVID-19_Related_Gene_Sets", destfile = "covid-19.gmt")
ln = readLines("covid-19.gmt")
ln = strsplit(ln, "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])
gs[1:2]

## -----------------------------------------------------------------------------
df = data.frame(
    geneset = rep(names(gs), times = sapply(gs, length)),
    gene = unlist(gs)
)
head(df)

## -----------------------------------------------------------------------------
df = list_to_dataframe(gs)
head(df)

## ----echo = FALSE-------------------------------------------------------------
invisible(file.remove("covid-19.gmt"))

