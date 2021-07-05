.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Accessory functions for RPP. Enjoy!"
  )
}

.onLoad <- function(libname, pkgname) {
  
  ass = function(name, x) {
    # double parent.env to get out of function body
    assign(name, x, envir = parent.env(parent.env(environment())))
  }

  ## Colorblind friendly palettes

  # The palette with grey:
  ass('cbPalette', c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73",
  "#D55E00", "#F0E442", "#0072B2"))
  
  # The palette with black:
  ass('cbbPalette', c("#000000", "#E69F00", "#56B4E9", "#CC79A7", "#009E73",
  "#D55E00", "#F0E442", "#0072B2"))

  ass('drugCols', 
    c(RifMono = '#a40000', PmbMono = '#e69e00', Cmb = '#000000',
      RifWt   = '#ef2929', PmbWt = '#fcaf3e', CmbWt = '#888a85',
      LoeweNullWt = '#babdb6', LoeweNullMut = '#2e3436'
     )
  )
  
  ass('myBreaks', c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28))
  
  ass('myCols', c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
      "#B15928"))

}
