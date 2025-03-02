## -----------------------------------------------------------------------------
library(Orthology.eg.db)

## -----------------------------------------------------------------------------
Orthology.eg.db

## -----------------------------------------------------------------------------
cl = columns(Orthology.eg.db)
length(cl)
head(cl)

## -----------------------------------------------------------------------------
kt = keytypes(Orthology.eg.db)
identical(cl, kt)

## -----------------------------------------------------------------------------
kt[grep("Ailuropoda", kt, ignore.case = TRUE)]

## -----------------------------------------------------------------------------
keys = keys(Orthology.eg.db, keytype = "Homo.sapiens")
head(keys)
length(keys)

## -----------------------------------------------------------------------------
map_df = select(Orthology.eg.db, keys, 
    columns = "Ailuropoda.melanoleuca", keytype = "Homo.sapiens")
head(map_df)

## -----------------------------------------------------------------------------
class(map_df[, 1])
class(map_df[, 2])

## -----------------------------------------------------------------------------
map_df[, 1] = as.character(map_df[, 1])
map_df[, 2] = as.character(map_df[, 2])

## -----------------------------------------------------------------------------
map_vec = structure(map_df[, 2], names = map_df[, 1])
# or we can do it in two lines
map_vec = map_df[, 2]
names(map_vec) = map_df[, 1]

head(map_vec)

## -----------------------------------------------------------------------------
library(GSEAtraining)
gs_human = get_msigdb(version = "2023.2.Hs", collection = "h.all")

## -----------------------------------------------------------------------------
gs_panda = lapply(gs_human, function(x) {
    x2 = map_vec[x]
    unname(x2[!is.na(x2)])
})

## -----------------------------------------------------------------------------
gs_panda = gs_panda[sapply(gs_panda, length) > 0]
gs_panda[1:2]

## -----------------------------------------------------------------------------
kt[grep("Tursiops", kt, ignore.case = TRUE)]

## -----------------------------------------------------------------------------
map_df = select(Orthology.eg.db, keys, 
    columns = "Tursiops.truncatus", keytype = "Homo.sapiens")
map_df[, 1] = as.character(map_df[, 1])
map_df[, 2] = as.character(map_df[, 2])

map_vec = structure(map_df[, 2], names = map_df[, 1])

gs_dolphin = lapply(gs_human, function(x) {
    x2 = map_vec[x]
    unname(x2[!is.na(x2)])
})
gs_dolphin = gs_dolphin[sapply(gs_dolphin, length) > 0]
gs_dolphin[1:2]

