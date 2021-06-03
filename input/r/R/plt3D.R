#' Plot drug-drug interaction plots in 3D
#'
#' Mostly a convenience function to have the most common settings/appearance for
#' manuscript: viewangle, consistent colors, range of z-axis etc.
#' @export
plt3D = function(mod, zlim=c(-5.5, 5), xlim=NULL, ylim=NULL){

  do.call(pltRS, list(mod, radius = 0.13, main = '', legend = F,
    zlab = '', xlab='', ylab = '', zlim = zlim, xlim=xlim, ylim=ylim))
 
  # View angle 
  rgl_graph_orientation = dget("input/dat/rgl_orientation")
  rgl.viewpoint(userMatrix = rgl_graph_orientation)
  
  # Add line segments to feature the difference
  interleave = function(v1, v2)  as.vector(rbind(v1,v2))
  tmp = mod$offAxisTable
  my_cols = ifelse(tmp$effect - tmp$predicted>0, "#2166AC", "#B2182B")
  log10T = function(z) log10(z + 0.5 * min(z[z != 0]))
  segments3d(interleave(log10T(tmp$d1), log10T(tmp$d1)),
             interleave(log10T(tmp$d2), log10T(tmp$d2)),
             interleave(tmp$effect,  tmp$predicted),
             lwd=1.5, col=interleave(my_cols, my_cols))

  aspect3d(x=1, y=1, z=1)
}
