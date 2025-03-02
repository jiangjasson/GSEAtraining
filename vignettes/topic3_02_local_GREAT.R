## -----------------------------------------------------------------------------
library(rGREAT)
set.seed(123)
gr = randomRegions(nr = 1000, genome = "hg19")

## ----echo = FALSE-------------------------------------------------------------
great_opt$verbose = FALSE

## -----------------------------------------------------------------------------
res = great(gr, "MSigDB:H", "hg19")

## -----------------------------------------------------------------------------
tb = getEnrichmentTable(res)
plotVolcano(res)
plotRegionGeneAssociations(res)
getRegionGeneAssociations(res)
plotRegionGeneAssociations(res, term_id = "HALLMARK_APOPTOSIS")
# shinyReport(res)

## -----------------------------------------------------------------------------
gs = read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"), 
    from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
gs[1:2]
great(gr, gs, "hg19")

## -----------------------------------------------------------------------------
df = read.table(url("https://rest.kegg.jp/link/pathway/hsa"), sep = "\t")
df[, 1] = gsub("hsa:", "", df[, 1])
df[, 2] = gsub("path:",  "", df[, 2])
gs_kegg = split(df[, 1], df[, 2])
gs_kegg[1:2]
great(gr, gs_kegg, "hg19")

## -----------------------------------------------------------------------------
library(reactome.db)
gs_reactome = as.list(reactomePATHID2EXTID)
great(gr, gs_reactome, "hg19")

## -----------------------------------------------------------------------------
tb = rGREAT:::BIOC_ANNO_PKGS
knitr::kable(tb[!duplicated(tb$genome_version_in_txdb), ])

## -----------------------------------------------------------------------------
library(msigdbr)
msigdbr_species()

## -----------------------------------------------------------------------------
h_gene_sets = msigdbr(species = "chimpanzee", category = "H")
head(h_gene_sets)

## -----------------------------------------------------------------------------
h_gene_sets = split(h_gene_sets$entrez_gene, h_gene_sets$gs_name)
h_gene_sets = lapply(h_gene_sets, as.character)  # just to make sure gene IDs are all in character.
h_gene_sets = lapply(h_gene_sets, function(x) unique(x[!is.na(x)])) # remove NA and duplicated genes
h_gene_sets[1:2]

## -----------------------------------------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()
query(ah, "TxDb")

