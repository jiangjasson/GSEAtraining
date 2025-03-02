## -----------------------------------------------------------------------------
library(UniProtKeywords)

## -----------------------------------------------------------------------------
UniProtKeywords

## -----------------------------------------------------------------------------
gl = load_keyword_genesets(9606)
gl[3:4]  # because gl[1:2] has a very long output, here we print gl[3:4]

## ----eval = FALSE-------------------------------------------------------------
# load_keyword_genesets("human")
# load_keyword_genesets("Homo sapiens")

## -----------------------------------------------------------------------------
tb = load_keyword_genesets(9606, as_table = TRUE)
head(tb)

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot(table(sapply(gl, length)), log = "x", 
    xlab = "Size of keyword genesets",
    ylab = "Number of keywords"
)

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot(table(sapply(gregexpr(" |-|/", names(gl)), length)), 
    xlab = "Number of words in keywords",
    ylab = "Number of keywords"
)

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot(table(nchar(names(gl))), 
    xlab = "Number of characters in keywords",
    ylab = "Number of keywords"
)

