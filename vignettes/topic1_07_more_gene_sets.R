## ----eval = FALSE-------------------------------------------------------------
# read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
# read.table(url("https://rest.kegg.jp/list/pathway/ptr"), sep = "\t")
# read.table(url("https://rest.kegg.jp/list/pathway/pps"), sep = "\t")
# read.table(url("https://rest.kegg.jp/list/pathway/ggo"), sep = "\t")
# read.table(url("https://rest.kegg.jp/list/pathway/pon"), sep = "\t")
# ...

## -----------------------------------------------------------------------------
library(BioMartGOGeneSets)
lt = getBioMartGOGeneSets("mmusculus_gene_ensembl")
length(lt)
lt[1]

## ----eval = FALSE-------------------------------------------------------------
# lt = getBioMartGOGeneSets("mouse")

## -----------------------------------------------------------------------------
tb = getBioMartGOGeneSets("mmusculus_gene_ensembl", as_table = TRUE)
head(tb)

## ----eval = FALSE-------------------------------------------------------------
# getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "BP") # the default one
# getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "CC")
# getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "MF")

## -----------------------------------------------------------------------------
lt = getBioMartGOGeneSets("mmusculus_gene_ensembl", gene_id_type = "entrez_gene")
lt[1]

lt = getBioMartGOGeneSets("mmusculus_gene_ensembl", gene_id_type = "gene_symbol")
lt[1]

## ----echo = FALSE-------------------------------------------------------------
tb = read.table(textConnection(
"org.Hs.eg.db    Human   org.Ss.eg.db    Pig
org.Mm.eg.db    Mouse   org.Gg.eg.db    Chicken
org.Rn.eg.db    Rat org.Mmu.eg.db   Rhesus_monkey
org.Dm.eg.db    Fruit_fly   org.Cf.eg.db    Canine
org.At.tair.db  Arabidopsis org.EcK12.eg.db E_coli_strain_K12
org.Sc.sgd.db   Yeast   org.Xl.eg.db    African_clawed_frog
org.Dr.eg.db    Zebrafish   org.Ag.eg.db    Malaria_mosquito
org.Ce.eg.db    Nematode    org.Pt.eg.db    Chimpanzee
org.Bt.eg.db    Bovine  org.EcSakai.eg.db   E_coli_strain_Sakai
"))

## ----echo = FALSE-------------------------------------------------------------
tb[, 1] = paste0("`", tb[, 1], "`")
tb[, 3] = paste0("`", tb[, 3], "`")
tb[, 2] = gsub("_", " ", tb[, 2])
tb[, 4] = gsub("_", " ", tb[, 4])
knitr::kable(tb, col.names = c("Package", "Organism", "Package", "Organism"))

## -----------------------------------------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()

## -----------------------------------------------------------------------------
query(ah, c("cat", "OrgDb"))

## -----------------------------------------------------------------------------
query(ah, c("Felis catus", "OrgDb"))

## -----------------------------------------------------------------------------
# It is annoying that ID changes between different package versions
org_db = ah[["AH117537"]]  # using `[[` downloads the dataset

## -----------------------------------------------------------------------------
org_db
columns(org_db)

## -----------------------------------------------------------------------------
all_genes = keys(org_db, keytype = "ENTREZID")
tb = select(org_db, keys = all_genes, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL"))
head(tb)

## -----------------------------------------------------------------------------
tb = tb[!is.na(tb$GOALL), ]
tb = unique(tb)

## -----------------------------------------------------------------------------
library(GSEAtraining)
lt = get_GO_gene_sets_from_orgdb(org_db, "BP")
lt[1]

## -----------------------------------------------------------------------------
gs1 = getBioMartGOGeneSets("ttruncatus_gene_ensembl", "BP", gene_id_type = "entrez_id")

## -----------------------------------------------------------------------------
qu = query(ah, c("Tursiops truncatus", "OrgDb"))
orgdb = ah[[ qu$ah_id ]]
gs2 = get_GO_gene_sets_from_orgdb(orgdb, "BP")

## -----------------------------------------------------------------------------
cn = intersect(names(gs1), names(gs2))
n1 = sapply(gs1, length)
n2 = sapply(gs2, length)

plot(n1[cn], n2[cn], xlab = "getBioMartGOGeneSets", ylab = "AnnotationHub")

