
#' GSEA plot
#' 
#' @param s A numeric vector of gene scores.
#' @param which Which gene sets to draw. Multiple gene sets are allowed.
#' @param gs A list of gene sets. Genes should have the same gene ID type as in `s`.
#' @param col Colors.
#' @param panel_height Relative height of the three panels in the plot.
#' 
#' @export
#' @import ComplexHeatmap
#' @import grid
gsea_plot = function(s, which, gs, col = seq_along(which), panel_height = c(6, 1, 3)) {

	s = sort(s, decreasing = TRUE)

	get_es = function(s, l_set) {
	    
	    s_set = abs(s)
	    s_set[!l_set] = 0

	    l_other = !l_set

	    f1 = cumsum(s_set)/sum(s_set)
	    f2 = cumsum(l_other)/sum(l_other)

	    f1 - f2
	}

	llt = lapply(which, function(i) {
		names(s) %fin% gs[[i]]
	})

	esl = lapply(llt, function(l) {
		get_es(s, l)
	})
    
    rg = range(unlist(lapply(esl, range)))
    rg = rg + c(-1, 1)*diff(rg)*0.05
    rg[1] = min(rg[1], 0)
    rg[2] = max(rg[2], 0)
    n = length(s)
    k = length(esl)

    panel_height = panel_height/sum(panel_height)

    lgd = Legend(title = "", at = names(gs[which]), type = "lines", legend_gp = gpar(col = col), direction = "horizontal",
    	background = NA, ncol = 3, labels_gp = gpar(fontsize = 8), by_row = TRUE)
    lgd_height = ComplexHeatmap:::height(lgd)
    h = unit(1, "npc") - unit(1 + 0.5 + 0.5 + 1.5, "cm") - lgd_height - unit(2, "mm")
    
    grid.newpage()
    pushViewport(viewport(width = unit(1, "npc") - unit(2, "cm"), height = h*panel_height[1], 
    	x = unit(1.5, "cm"), y = unit(1, "npc") - unit(1, "cm"), 
    	just = c("left", "top"), xscale = c(0, n+1), yscale = rg))
    grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "grey", lty = 2))
    for(i in seq_along(esl)) {
    	grid.lines(1:n, esl[[i]], gp = gpar(col = col[i]), default.units = "native")
    	j = which.max(abs(esl[[i]]))
    	grid.lines(c(j, j), c(esl[[i]][j], 0), default.units = "native", gp = gpar(lty = 2, col = col[i]))
    }
    grid.rect()
    grid.yaxis(at = grid.pretty(rg), gp = gpar(fontsize = 8))
    grid.text("Enrichment scores", x = unit(-1.2, "cm"), y = unit(0.5, "npc"), rot = 90)
    popViewport()

    pushViewport(viewport(width = unit(1, "npc") - unit(2, "cm"), height = h*panel_height[2], 
    	x = unit(1.5, "cm"), y = unit(1, "npc") - unit(1, "cm") - h*panel_height[1] - unit(0.5, "cm"), 
    	just = c("left", "top"), xscale = c(0, n+1), yscale = c(0, k)))
    for(i in seq_along(esl)) {
    	grid.segments(which(llt[[i]]), i-0.8, which(llt[[i]]), i-0.2, default.units = "native", gp = gpar(col = col[i]))
    }
    grid.rect()
    grid.text("Gene sets", x = unit(-1.2, "cm"), y = unit(0.5, "npc"), rot = 90)
    popViewport()

    pushViewport(viewport(width = unit(1, "npc") - unit(2, "cm"), height = h*panel_height[3], 
    	x = unit(1.5, "cm"), y = unit(1, "npc") - unit(1, "cm") - h*panel_height[1] - unit(0.5, "cm") - h*panel_height[2] - unit(0.5, "cm"),  
    	just = c("left", "top"), xscale = c(0, n+1), yscale = range(s)))
    grid.polygon(c(1, 1:n, n), c(0, s, 0), default.units = "native", gp = gpar(fill = "grey", col = "grey"))
    grid.rect(gp = gpar(fill = NA))
    grid.yaxis(at = grid.pretty(range(s), n = 3), gp = gpar(fontsize = 8))
    grid.text("Gene scores", x = unit(-1.2, "cm"), y = unit(0.5, "npc"), rot = 90)
    grid.xaxis(gp = gpar(fontsize = 8))
    grid.text("Ranked gene list", y = unit(-1, "cm"))
    popViewport()
    
    pushViewport(viewport(width = unit(1, "npc") - unit(2, "cm"), height = lgd_height, 
    	x = unit(1.5, "cm"), y = unit(1, "npc") - unit(1, "cm") - h*panel_height[1] - unit(0.5, "cm") - h*panel_height[2] - unit(0.5, "cm") - h*panel_height[3] - unit(1.5, "cm"),  
    	just = c("left", "top"), ))
    grid.draw(lgd)
    popViewport()

}
