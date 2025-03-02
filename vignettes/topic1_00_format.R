## -----------------------------------------------------------------------------
lt = list(
    geneset1 = c("gene1", "gene2", "gene3"),
    geneset2 = c("gene2", "gene4"),
    geneset3 = c("gene1", "gene3", "gene5", "gene6")
)
lt

## -----------------------------------------------------------------------------
df = data.frame(
    geneset = c(rep("geneset1", 3), rep("geneset2", 2), rep("geneset3", 4)),
    gene = c("gene1", "gene2", "gene3", "gene2", "gene4", "gene1", "gene3", "gene5", "gene6")
)
df

## -----------------------------------------------------------------------------
df[, 2:1]

## -----------------------------------------------------------------------------
data.frame(
    geneset = rep(names(lt), times = sapply(lt, length)),
    gene = unlist(lt)
)

## -----------------------------------------------------------------------------
split(df$gene, df$geneset)

## -----------------------------------------------------------------------------
library(GSEAtraining)
list_to_dataframe(lt)
dataframe_to_list(df)

## -----------------------------------------------------------------------------
m = matrix(0, nrow = 3, ncol = 6)
rownames(m) = unique(df$geneset)
colnames(m) = unique(df$gene)

for(i in seq_len(nrow(df))) {
    m[df[i, 1], df[i, 2]] = 1
}
m

