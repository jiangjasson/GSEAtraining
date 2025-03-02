## -----------------------------------------------------------------------------
library(BioCartaImage)
head(PATHWAY2ENTREZ)

## -----------------------------------------------------------------------------
pathway_names = sapply(BIOCARTA_PATHWAYS, function(x) x$name)
head(pathway_names)

## -----------------------------------------------------------------------------
library(grid)
grid.newpage()
grid.biocarta("h_RELAPathway", color = c("1387" = "yellow"))

## -----------------------------------------------------------------------------
grob = biocartaGrob("h_RELAPathway")
grob2 = mark_gene(grob, "1387", function(x, y) {
    pos = pos_by_polygon(x, y)
    pushViewport(viewport(x = pos[1] - 10, y = pos[2], 
        width = unit(4, "cm"), height = unit(4, "cm"), 
        default.units = "native", just = "right"))
    grid.rect(gp = gpar(fill = "red"))
    grid.text("add whatever\nyou want here")
    popViewport()
}, capture = TRUE)
grid.draw(grob2)

