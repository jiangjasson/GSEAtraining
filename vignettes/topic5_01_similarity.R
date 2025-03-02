## -----------------------------------------------------------------------------
library(GSEAtraining)
lt = readRDS(system.file("extdata", "ora.rds", package = "GSEAtraining"))
diff_gene = lt$diff_gene
diff_gene = convert_to_entrez_id(diff_gene)

## -----------------------------------------------------------------------------
df = read.table(url("https://rest.kegg.jp/link/pathway/hsa"), sep = "\t")
df[, 1] = gsub("hsa:", "", df[, 1])
df[, 2] = gsub("path:",  "", df[, 2])
head(df)

## -----------------------------------------------------------------------------
pathway_names = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
pathway_names = structure(gsub(" - Homo .*$", "", pathway_names[, 2]), names = pathway_names[, 1])

## -----------------------------------------------------------------------------
library(clusterProfiler)
tb = enricher(diff_gene, TERM2GENE = df[, 2:1])
sig_pathways = as.data.frame(tb)$ID

gs_pathway = split(df[, 1], df[, 2])
gs_pathway = gs_pathway[sig_pathways]

## -----------------------------------------------------------------------------
m1 = term_similarity(gs_pathway, method = "jaccard")
m2 = term_similarity(gs_pathway, method = "overlap")
m3 = term_similarity(gs_pathway, method = "kappa")

## ----fig.width = 10, fig.height = 7-------------------------------------------
diag(m1) = NA
diag(m2) = NA
diag(m3) = NA

library(ComplexHeatmap)
Heatmap(m1, name = "jaccard") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))
Heatmap(m2, name = "overlap") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))
Heatmap(m3, name = "kappa") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))

## ----fig.width = 10, fig.height = 7-------------------------------------------
m1 = term_similarity(gs_pathway, method = "jaccard", all = diff_gene)
m2 = term_similarity(gs_pathway, method = "overlap", all = diff_gene)
m3 = term_similarity(gs_pathway, method = "kappa", all = diff_gene)

diag(m1) = NA
diag(m2) = NA
diag(m3) = NA

Heatmap(m1, name = "jaccard") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))
Heatmap(m2, name = "overlap") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))
Heatmap(m3, name = "kappa") + 
    rowAnnotation(pathway_name = anno_text(pathway_names[rownames(m1)]))

## -----------------------------------------------------------------------------
library(org.Hs.eg.db)
tb = enrichGO(gene = diff_gene, ont = "BP", OrgDb = org.Hs.eg.db)

sig_go = as.data.frame(tb)$ID

## -----------------------------------------------------------------------------
gs_bp = get_GO_gene_sets_from_orgdb(org.Hs.eg.db)
gs_bp = gs_bp[sig_go]

## -----------------------------------------------------------------------------
m1 = term_similarity(gs_bp, method = "jaccard")
m2 = term_similarity(gs_bp, method = "overlap")
m3 = term_similarity(gs_bp, method = "kappa")

diag(m1) = NA
diag(m2) = NA
diag(m3) = NA

Heatmap(m1, name = "jaccard", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(m2, name = "overlap", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(m3, name = "kappa", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)

## -----------------------------------------------------------------------------
m1 = term_similarity(gs_bp, method = "jaccard", all = diff_gene)
m2 = term_similarity(gs_bp, method = "overlap", all = diff_gene)
m3 = term_similarity(gs_bp, method = "kappa", all = diff_gene)

diag(m1) = NA
diag(m2) = NA
diag(m3) = NA

Heatmap(m1, name = "jaccard", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(m2, name = "overlap", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(m3, name = "kappa", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)

## -----------------------------------------------------------------------------
library(simona)
simona_opt$verbose = FALSE

## -----------------------------------------------------------------------------
dag = create_ontology_DAG_from_GO_db("BP", org_db = org.Hs.eg.db)

## -----------------------------------------------------------------------------
m = term_sim(dag, sig_go, method = "Sim_Lin_1998")
diag(m) = NA
Heatmap(m, name = "Sim_Lin_1998", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)

## -----------------------------------------------------------------------------
m = term_sim(dag, sig_go, method = "Sim_WP_1994")
diag(m) = NA
Heatmap(m, name = "Sim_WP_1994", show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE)

## ----fig.width = 10, fig.height = 7-------------------------------------------
dag_circular_viz(dag, highlight = sig_go)

