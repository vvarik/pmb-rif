#' Plot drug-drug interaction plots in 3D
#'
#' Mostly a convenience function to have the most common settings/appearance for
#' manuscript: viewangle, consistent colors, range of z-axis etc.
#' @export
plt3D = function(mod, zlim=c(-4.5, 4)){

  do.call(pltRS, list(mod, radius = 0.1, main = '', legend = F,
    zlab = '', xlab='', ylab = '', zlim = zlim))
 
  # View angle 
  rgl_graph_orientation = dget("input/dat/rgl_orientation")
  rgl.viewpoint(userMatrix = rgl_graph_orientation)
  
  # Add line segments to feature the difference
  interleave = function(v1, v2)  as.vector(rbind(v1,v2))
  tmp = mod$offAxisTable
  my_cols = ifelse(tmp$effect - tmp$predicted>0, "#2166AC", "#B2182B")
  segments3d(interleave(log10(tmp$d1), log10(tmp$d1)),
             interleave(log10(tmp$d2), log10(tmp$d2)),
             interleave(tmp$effect,  tmp$predicted),
             lwd=2, col=interleave(my_cols, my_cols))
}
