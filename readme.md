


Install all suggested dependencies:

```r
setRepositories(ind = 1:4) # to include both CRAN and bioc repos

install.packages(c("knitr", "rmarkdown", "GO.db", "org.Hs.eg.db", "KEGGREST", 
	"clusterProfiler", "msigdbr", "reactome.db", "BioCartaImage", "UniProtKeywords", 
	"BioMartGOGeneSets", "AnnotationHub", "BiocHubsShiny", "microbenchmark", "ReactomePA", 
	"DOSE", "org.Ss.eg.db", "CePa", "eulerr", "rGREAT", "goseq", "GSVA", "simplifyEnrichment", 
	"simona", "enrichplot", "ggplot2", "ComplexHeatmap", "circlize", "genefilter", 
	"cola", "proxyC", "ggupset", "ggridges"))
```

**Orthology.eg.db** in the current bioc version (3.20) seems to have a problem. We use a lower version:

```r
install.packages("https://www.bioconductor.org/packages/3.17/data/annotation/src/contrib/Orthology.eg.db_3.17.0.tar.gz", 
	repo = NULL, type = "source")
```

Then

```r
install.packages("https://jokergoo.github.io/GSEAtraining/GSEAtraining_3.20.0.tar.gz", 
	repo = NULL, type = "source")
```

You may need to update the **rGREAT** package because the API link from NCBI was changed:

```r
library(devtools)
install_github("jokergoo/rGREAT")
```

The practice materials are also available at https://jokergoo.github.io/GSEAtraining/.

**Do not use it for commercial purpose.**

