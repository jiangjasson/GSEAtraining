## -----------------------------------------------------------------------------
ora = function(genes, gene_sets, universe = NULL) {

    if(is.null(universe)) {
        universe = unique(unlist(gene_sets))
    } else {
        universe = unique(universe)
    }
    # restrict in universe
    genes = intersect(genes, universe)
    gene_sets = lapply(gene_sets, function(x) intersect(x, universe))

    n_universe = length(universe)
    n_genes = length(genes)
    
    x = sapply(gene_sets, function(x) length(intersect(x, genes)))
    m = sapply(gene_sets, length)
    n = n_universe - m
    k = n_genes
    
    p = phyper(x - 1, m, n, k, lower.tail = FALSE)

    df = data.frame(
        gene_set = names(gene_sets),
        hits = x,
        n_genes = k,
        n_gs = m,
        n_total = n_universe,
        log2fe = log2(x*n_universe/k/m),
        p_value = p
    )

    df$p_adjust = p.adjust(df$p_value, "BH")
    rownames(df) = df$gene_set
    df[order(df$p_adjust, df$p_value), ,drop = FALSE]
}

## -----------------------------------------------------------------------------
library(org.Hs.eg.db)
library(GSEAtraining)
gs = get_msigdb(version = "2023.2.Hs", collection = "h.all")
genes = random_genes(org.Hs.eg.db, 1000, "ENTREZID")

df = ora(genes, gs)
head(df)

## -----------------------------------------------------------------------------
de = readRDS(system.file("extdata", "de.rds", package = "GSEAtraining"))
de = de[, c("symbol", "p_value")]
de = de[!is.na(de$p_value), ]
head(de)

## -----------------------------------------------------------------------------
de_genes_1 = de$symbol[de$p_value < 0.05]
de_genes_2 = de$symbol[de$p_value < 0.01]
de_genes_3 = de$symbol[de$p_value < 0.001]

length(de_genes_1)
length(de_genes_2)
length(de_genes_3)

## -----------------------------------------------------------------------------
gs = get_GO_gene_sets_from_orgdb(org.Hs.eg.db, "BP", gene_id_type = "SYMBOL")
gs[1]

## -----------------------------------------------------------------------------
tb1 = ora(de_genes_1, gs)
tb2 = ora(de_genes_2, gs)
tb3 = ora(de_genes_3, gs)

## -----------------------------------------------------------------------------
sum(tb1$p_adjust < 0.05)
sum(tb2$p_adjust < 0.05)
sum(tb3$p_adjust < 0.05)

## -----------------------------------------------------------------------------
library(eulerr)
plot(euler(list("DE_cutoff_0.05"  = tb1$gene_set[tb1$p_adjust < 0.05],
                "DE_cutoff_0.01"  = tb2$gene_set[tb2$p_adjust < 0.05],
                "DE_cutoff_0.001" = tb3$gene_set[tb3$p_adjust < 0.05])), 
     quantities = TRUE)

## -----------------------------------------------------------------------------
cn = intersect(intersect(rownames(tb1), rownames(tb2)), rownames(tb3))

## ----fig.width = 12, fig.height = 4-------------------------------------------
par(mfrow = c(1, 3))
plot(tb1[cn, "log2fe"], tb2[cn, "log2fe"], pch = 16, cex = 0.5, col = "#00000080",
    xlim = c(-5, 5), ylim = c(-5, 5), main = "compare log2 fold enrichment from ORA",
    xlab = "DE_cutoff_0.05", ylab = "DE_cutoff_0.01")
abline(a = 0, b = 1, lty = 2, col = "red")

plot(tb1[cn, "log2fe"], tb3[cn, "log2fe"], pch = 16, cex = 0.5, col = "#00000080",
    xlim = c(-5, 5), ylim = c(-5, 5), main = "compare log2 fold enrichment from ORA",
    xlab = "DE_cutoff_0.05", ylab = "DE_cutoff_0.001")
abline(a = 0, b = 1, lty = 2, col = "red")

plot(tb2[cn, "log2fe"], tb3[cn, "log2fe"], pch = 16, cex = 0.5, col = "#00000080",
    xlim = c(-5, 5), ylim = c(-5, 5), main = "compare log2 fold enrichment from ORA",
    xlab = "DE_cutoff_0.01", ylab = "DE_cutoff_0.001")
abline(a = 0, b = 1, lty = 2, col = "red")

