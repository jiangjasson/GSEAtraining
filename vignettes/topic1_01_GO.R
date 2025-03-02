## ----message = FALSE----------------------------------------------------------
library(GO.db)

## -----------------------------------------------------------------------------
GO.db

## -----------------------------------------------------------------------------
select(GO.db, keys = c("GO:0000001", "GO:0000002"), columns = c("ONTOLOGY", "TERM"))

## -----------------------------------------------------------------------------
columns(GO.db)

## -----------------------------------------------------------------------------
lt = as.list(GOBPCHILDREN)
head(lt)

## -----------------------------------------------------------------------------
tb = toTable(GOBPCHILDREN)
head(tb)

## -----------------------------------------------------------------------------
colnames(tb)[1:2] = c("child", "parent")

## -----------------------------------------------------------------------------
tb = toTable(GOBPCHILDREN)
table(tb[, 3])
barplot(table(tb[, 3]), main = "BP")

tb = toTable(GOMFCHILDREN)
table(tb[, 3])
barplot(table(tb[, 3]), main = "MF")

tb = toTable(GOCCCHILDREN)
table(tb[, 3])
barplot(table(tb[, 3]), main = "CC")

## -----------------------------------------------------------------------------
lt = as.list(GOBPCHILDREN)
tb = table(sapply(lt, length))
x = as.numeric(names(tb))
y = as.vector(tb)

library(ggplot2)
ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point() +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(x = "Number of child terms", y = "Number of GO terms") + ggtitle("GOBPCHILDREN")

## -----------------------------------------------------------------------------
lt = as.list(GOBPPARENTS)
lt = lt[names(lt) != "GO:0008150"]
tb = table(sapply(lt, length))
x = as.numeric(names(tb))
y = as.vector(tb)

ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point() +
    scale_y_continuous(trans='log10') +
    labs(x = "Number of parent terms", y = "Number of GO terms") + ggtitle("GOBOPARENTS")

## -----------------------------------------------------------------------------
lt = as.list(GOBPOFFSPRING)
lt[1:2]
tb = toTable(GOBPOFFSPRING)
head(tb)

## -----------------------------------------------------------------------------
GOBPCHILDREN[["GO:0000002"]]
GOBPOFFSPRING[["GO:0000002"]]

## -----------------------------------------------------------------------------
GOTERM

## -----------------------------------------------------------------------------
lt = as.list(GOTERM)
tb = toTable(GOTERM)

## -----------------------------------------------------------------------------
head(Term(GOTERM))
head(Definition(GOTERM))
head(Ontology(GOTERM))

## -----------------------------------------------------------------------------
GOTERM[c("GO:0000001", "GO:0000002", "GO:0000006")]

## -----------------------------------------------------------------------------
Term(GOTERM[c("GO:0000001", "GO:0000002", "GO:0000006")])
Definition(GOTERM[c("GO:0000001", "GO:0000002", "GO:0000006")])
Ontology(GOTERM[c("GO:0000001", "GO:0000002", "GO:0000006")])

## -----------------------------------------------------------------------------
library(org.Hs.eg.db)

## -----------------------------------------------------------------------------
org.Hs.eg.db

## -----------------------------------------------------------------------------
lt = as.list(org.Hs.egGO2ALLEGS)
lt[3:4]

## -----------------------------------------------------------------------------
lt2 = as.list(org.Hs.egGO2EG)
lt2[3:4]

## -----------------------------------------------------------------------------
tb = toTable(org.Hs.egGO2ALLEGS)
head(tb)

## -----------------------------------------------------------------------------
barplot(sort(table(tb$Evidence)))

## -----------------------------------------------------------------------------
org.Hs.egGO2ALLEGS[["GO:0000002"]]

## -----------------------------------------------------------------------------
n1 = sapply(lt, length)
lt = lapply(lt, unique)
n2 = sapply(lt, length)

plot(n1, n2, xlim = c(0, max(n1)), ylim = c(0, max(n1)),
    xlab = "with duplicated", ylab = "without duplicated",
    main = "number of genes in GO BP gene sets")

tb = tb[, 1:2]
tb = unique(tb)

## -----------------------------------------------------------------------------
tb = toTable(org.Hs.egGO2ALLEGS)
tb = unique(tb[tb$Ontology == "BP", 1:2])
t1 = table(table(tb$go_id))
x1 = as.numeric(names(t1))
y1 = as.vector(t1)
ggplot(data.frame(x = x1, y = y1), aes(x = x, y = y)) +
    geom_point() +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(x = "Number of annotated genes", y = "Number of GO terms") + ggtitle("GOBP")

## -----------------------------------------------------------------------------
t2 = table(table(tb$gene_id))
x2 = as.numeric(names(t2))
y2 = as.vector(t2)
ggplot(data.frame(x = x2, y = y2), aes(x = x, y = y)) +
    geom_point() +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(x = "Number of gene sets", y = "Number of genes") + ggtitle("GOBP")

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
all_genes = keys(org.Hs.eg.db, keytype = "ENTREZID")
tb = select(org.Hs.eg.db, keys = all_genes, keytype = "ENTREZID", 
    columns = c("GOALL", "ONTOLOGYALL"))
head(tb)

## -----------------------------------------------------------------------------
tb = tb[tb$ONTOLOGYALL == "BP", 1:2]
tb = unique(tb)
head(tb)

## -----------------------------------------------------------------------------
lt = as.list(org.Hs.egGO2ALLEGS)
n = length(lt)
nd = sum(sapply(lt, function(x) any(duplicated(x))))
nd/n

## -----------------------------------------------------------------------------
lt_genes = as.list(org.Hs.egGO2EG)
lt_terms = as.list(GOBPOFFSPRING)

## -----------------------------------------------------------------------------
for(nm in names(lt_terms)) {
    lt_terms[[nm]] = c(lt_terms[[nm]], nm)
}

## -----------------------------------------------------------------------------
lt_genes_manual = lapply(lt_terms, function(x) {
    unique(unlist(lt_genes[x]))
})

## -----------------------------------------------------------------------------
lt_genes_all = as.list(org.Hs.egGO2ALLEGS)
lt_genes_all = lapply(lt_genes_all, unique)

cn = intersect(names(lt_genes_manual), names(lt_genes_all))
plot(sapply(lt_genes_manual[cn], length), sapply(lt_genes_all[cn], length), log = "xy",
    xlab = "manual gene sets", ylab = "org.Hs.egGO2ALLEGS")

