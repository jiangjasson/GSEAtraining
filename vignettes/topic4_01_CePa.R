## -----------------------------------------------------------------------------
library(CePa)
data(PID.db)
data(gene.list)

## ----eval = FALSE-------------------------------------------------------------
# # around 5min
# res = cepa.ora.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)

## ----echo = FALSE-------------------------------------------------------------
res = readRDS(system.file("extdata", "cepa_res.rds", package = "GSEAtraining"))

## -----------------------------------------------------------------------------
plot(res, adj.method = "BH", only.sig = TRUE)

## -----------------------------------------------------------------------------
plot(res, "e2f_pathway", cen = "in.degree")
plot(res, "e2f_pathway", cen = "out.degree")

