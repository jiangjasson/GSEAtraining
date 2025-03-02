## -----------------------------------------------------------------------------
bed = read.table(system.file("extdata", "mm9_chipseq_peaks.bed", package = "GSEAtraining"))
head(bed)
dim(bed)

## ----message = FALSE----------------------------------------------------------
library(rGREAT)
job = submitGreatJob(bed, species = "mm9")
job

## -----------------------------------------------------------------------------
availableOntologies(job)

## -----------------------------------------------------------------------------
tbl = getEnrichmentTables(job)

## -----------------------------------------------------------------------------
names(tbl)
head(tbl[["GO Biological Process"]])

## -----------------------------------------------------------------------------
plotVolcano(job, ontology = "GO Biological Process")

## ----fig.width = 10-----------------------------------------------------------
plotRegionGeneAssociationGraphs(job)

## ----fig.width = 10-----------------------------------------------------------
plotRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")

## -----------------------------------------------------------------------------
getRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")

## ----eval = FALSE-------------------------------------------------------------
# shinyReport(job)

