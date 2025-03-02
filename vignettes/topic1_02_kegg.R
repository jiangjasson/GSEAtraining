## -----------------------------------------------------------------------------
df1 = read.table(url("https://rest.kegg.jp/link/pathway/hsa"), sep = "\t")
head(df1)

## ----eval = FALSE-------------------------------------------------------------
# temp = tempfile()
# download.file("https://rest.kegg.jp/link/pathway/hsa", destfile = temp)
# df1 = read.table(temp, sep = "\t")
# file.remove(temp)

## -----------------------------------------------------------------------------
df1[, 1] = gsub("hsa:", "", df1[, 1])
df1[, 2] = gsub("path:",  "", df1[, 2])
head(df1)

## -----------------------------------------------------------------------------
df2 = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
head(df2)
df2[, 2] = gsub(" - Homo.*$", "", df2[, 2])
head(df2)

## -----------------------------------------------------------------------------
library(KEGGREST)
gene2pathway = keggLink("pathway", "hsa")
head(gene2pathway)

## -----------------------------------------------------------------------------
df3 = data.frame(
    gene_id    = gsub("hsa:", "", names(gene2pathway)),
    pathway_id = gsub("path:", "", gene2pathway)
)
head(df3)

## -----------------------------------------------------------------------------
head(keggList("pathway", "hsa"))

## ----comment = ""-------------------------------------------------------------
genes = names(gene2pathway[gene2pathway == "path:hsa04911"])
diff_genes = sample(genes, 10)
settings = data.frame(
    genes = diff_genes, 
    colors = paste0(rep("orange", 10),",",rep("blue", 10)) # background,labels
)
write.table(settings, file = stdout(), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

## ----warnings = FALSE, message = FALSE----------------------------------------
lt = clusterProfiler::download_KEGG("hsa")
head(lt$KEGGPATHID2EXTID)
head(lt$KEGGPATHID2NAME)

## ----eval = FALSE-------------------------------------------------------------
# df = read.table(url("https://rest.kegg.jp/link/pathway/cfa"), sep = "\t")
# pathway2gene = KEGGREST::keggLink("pathway", "cfa")
# lt = clusterProfiler::download_KEGG("cfa")

## -----------------------------------------------------------------------------
df1 = read.table(url("https://rest.kegg.jp/link/pathway/hsa"), sep = "\t")
n_gene = table(df1[, 2])
hist(n_gene)

## -----------------------------------------------------------------------------
hist(n_gene, nc = 100)

## -----------------------------------------------------------------------------
which.max(n_gene)

## -----------------------------------------------------------------------------
df2 = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
df2[df2[, 1] == "hsa01100", ]

## -----------------------------------------------------------------------------
df2 = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
df2[, 2] = gsub(" - Homo .*$", "", df2[, 2])
pathway_id = df2[df2[, 2] == "Cell cycle", 1]
pathway_id

