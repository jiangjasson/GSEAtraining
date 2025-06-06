
This package aims to provide a comprehensive introduction on gene set enrichment analysis.

## Topics

**Topic 1**: gene set resources

- Topic 1-00: [Representation of gene sets in R](articles/topic1_00_format.html)
- Topic 1-01: [GO gene sets](articles/topic1_01_GO.html)
- Topic 1-02: [Get pathways from KEGG](articles/topic1_02_kegg.html)
- Topic 1-03: [Get pathways from MSigDB](articles/topic1_03_msigdb.html)
- Topic 1-04: [Get pathways from Reactome](articles/topic1_04_Reactome.html)
- Topic 1-05: [Get gene sets from UniProt Keywords](articles/topic1_05_UniProtKeywords.html)
- Topic 1-06: [Get GO/KEGG gene sets for other organisms](articles/topic1_06_more_gene_sets.html)
- Topic 1-07: [Generate gene sets for other organisms by mapping to orthologues](articles/topic1_07_ortholog.html)
- Topic 1-08: [Gene ID conversion](articles/topic1_08_gene_id.html)


**Topic 2**: ORA and GSEA

- Topic 2-01: [ORA on a single gene set](articles/topic2_01_ora_single_gene_set.html)
- Topic 2-02: [Implement ORA from scratch](articles/topic2_02_implement_ora.html)
- Topic 2-03: [Implement GSEA from scratch](articles/topic2_03_implement_gsea.html)
- Topic 2-04: [fgsea](articles/topic2_04_fgsea.html)
- Topic 2-05: [Compare ORA and GSEA](articles/topic2_05_compare_ORA_GSEA.html)

**Topic 3**: Application in genomics

- Topic 3-01: [Implement the GREAT method](articles/topic3_01_implement_GREAT.html)
- Topic 3-02: [Online GREAT analysis](articles/topic3_02_online_GREAT.html)
- Topic 3-03: [Local GREAT analysis](articles/topic3_03_local_GREAT.html)

**Topic 4**: Topology-based

- Topic 4-01: [Centrality-based pathway enrichment analysis](articles/topic4_01_centrality.html)

**Topic 5**: Single-sample GSEA

- Topic 5-01: [Single-sample GSEA](articles/topic5_01_single_sample.html)

**Topic 6**: Similarity, clustering and summarization

- Topic 6-01: [Similarities of terms](articles/topic6_01_similarity.html)
- Topic 6-02: [Cluster and summarize](articles/topic6_02_cluster_and_summarize.html)

**Topic 7**: Visualization

- Topic 7-01: [Visualization](articles/topic7_01_visualization.html)

**Topic 8**: GSEA framework

- Topic 8-01: [GSEA framework](articles/topic8_01_GSEA_framework.html)

## Install

```r
library(devtools)
install_github("jokergoo/GSEAtopics")
```

To install all necessary packages for running the vignettes:

```r
setRepositories(ind = 1:4)
install.packages(c("circlize", "reactome.db", "UniProtKeywords",
	"AnnotationHub", "Orthology.eg.db", "GSVA", "simona", "BiocManager", "CePa",
	"ggplot2", "rGREAT", "KEGGgraph", "proxyC", "sparseMatrixStats", "HilbertCurve", 
	"TxDb.Hsapiens.UCSC.hg19.knownGene", "cola", "cowplot", "eulerr", "golubEsets", 
	"hu6800.db", "microbenchmark", "preprocessCore", "simplifyEnrichment"))
```

There is also an simplified version

https://carpentries-incubator.github.io/bioc-rnaseq/07-gene-set-analysis.html


## License

Code is released under [the MIT licence](https://opensource.org/license/mit). Vignettes are under [the CC BY-NC-ND license](https://creativecommons.org/licenses/by-nc-nd/4.0/).

