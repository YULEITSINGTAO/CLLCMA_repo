# Colorblind-friendly cell_type_final palette (Paul Tol + Okabe-Ito, fixed order).
# Matches scientific-visualization skill: no rainbow interpolation; redundant-safe hues.

cell_type_palette <- function(types) {
  types <- sort(unique(as.character(types)))
  base <- c(
    # Paul Tol bright
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377",
    # Paul Tol vibrant
    "#EE7733", "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988",
    # Paul Tol muted
    "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255",
    "#44AA99", "#999933", "#AA4499",
    # Paul Tol medium contrast
    "#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA",
    # Paul Tol light (labeled fills)
    "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99",
    "#BBCC33", "#AAAA00",
    # Okabe-Ito
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#000000",
    # Paul Tol dark (sparingly)
    "#222255", "#225555", "#225522", "#666633", "#663333"
  )
  if (length(types) > length(base)) {
    stop("Need ", length(types), " colors but palette provides ", length(base))
  }
  stats::setNames(base[seq_along(types)], types)
}

label_text_color <- function(hex) {
  rgb <- grDevices::col2rgb(hex, alpha = FALSE) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum < 0.45, "#FFFFFF", "#1A1A1A")
}
